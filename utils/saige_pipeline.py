from hailtop import pipeline
from hailtop.pipeline.pipeline import *
from collections import Counter
from shlex import quote as shq


MKL_OFF = 'export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; '
SCRIPT_DIR = '/ukb_common/saige'
saige_pheno_types = {
    'continuous': 'quantitative',
    'biomarkers': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary',
    'icd_all': 'binary',
    'phecode': 'binary',
    'prescriptions': 'binary'
}


def create_sparse_grm(p: Pipeline, output_path: str, plink_file_root: str, docker_image: str,
                      relatedness_cutoff: str = '0.125', num_markers: int = 2000,
                      n_threads: int = 8, storage = '1500Mi'):
    in_bfile = p.read_input_group(
        **{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    create_sparse_grm_task: pipeline.pipeline.Task = p.new_task(name='create_sparse_grm')
    create_sparse_grm_task.cpu(n_threads).storage(storage).image(docker_image)
    create_sparse_grm_task.declare_resource_group(
        sparse_grm={ext: f'{{root}}{ext}' for ext in
                   (f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx',
                    f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt')})
    command = (f'Rscript /usr/local/bin/createSparseGRM.R '
               f'--plinkFile={in_bfile} '
               f'--nThreads={n_threads} '
               f'--outputPrefix={create_sparse_grm_task.sparse_grm}	'
               f'--numRandomMarkerforSparseKin={num_markers} '
               f'--relatednessCutoff={relatedness_cutoff}')
    create_sparse_grm_task.command(command)
    p.write_output(create_sparse_grm_task.sparse_grm, output_path)
    # Runtime: ~40 minutes on 8 cores
    return create_sparse_grm_task.sparse_grm


def extract_vcf_from_mt(p: Pipeline, output_root: str, docker_image: str, module: str = 'ukb_exomes',
                        gene: str = None, interval: str = None, groups=None, gene_map_ht_path: str = None,
                        set_missing_to_hom_ref: bool = False, callrate_filter: float = 0.0, adj: bool = True,
                        export_bgen: bool = True, input_dosage: bool = False, reference: str = 'GRCh38',
                        n_threads: int = 8, storage: str = '500Mi', additional_args: str = '', memory: str = ''):
    if groups is None:
        # groups = {'pLoF', 'missense|LC', 'pLoF|missense|LC', 'synonymous'}
        groups = {'pLoF', 'missense|LC', 'synonymous'}
    extract_task: pipeline.pipeline.Task = p.new_task(name='extract_vcf',
                                                      attributes={
                                                          'interval': interval
                                                      })
    extract_task.image(docker_image).cpu(n_threads).storage(storage)
    if export_bgen:
        extract_task.declare_resource_group(out={'bgen': '{root}.bgen',
                                                 'sample': '{root}.sample',
                                                 'bgen.bgi': '{root}.bgen.bgi'})
    else:
        extract_task.declare_resource_group(out={'vcf.gz': f'{{root}}.vcf.gz',
                                                 'vcf.gz.tbi': f'{{root}}.vcf.gz.tbi'})

    output_file = f'{extract_task.bgz}.bgz' if not export_bgen else extract_task.out
    command = f"""set -o pipefail; python3 {SCRIPT_DIR}/extract_vcf_from_mt.py
    --load_module {module}
    {"--additional_args " + additional_args if additional_args else ''}
    {"--gene " + gene if gene else ""}
    {"--interval " + interval if interval else ""}
    --groups "{','.join(groups)}"
    --reference {reference}
    {"--gene_map_ht_path " + gene_map_ht_path if gene_map_ht_path else ""} 
    {"--callrate_filter " + str(callrate_filter) if callrate_filter else ""} 
    {"--export_bgen" if export_bgen else ""}
    {"--input_bgen" if input_dosage else ""}
    {"" if set_missing_to_hom_ref else "--mean_impute_missing"}
    {"" if adj else "--no_adj"} 
    {"--group_output_file " + extract_task.group_file if gene_map_ht_path else ""}
    --output_file {output_file} | tee {extract_task.stdout}
    ;""".replace('\n', ' ')

    if export_bgen:
        command += f'\n/bgen_v1.1.4-Ubuntu16.04-x86_64/bgenix -g {extract_task.out.bgen} -index -clobber'
    else:
        command += f'\nmv {extract_task.bgz}.bgz {extract_task.out["vcf.gz"]}; tabix {extract_task.out["vcf.gz"]};'
    extract_task.command(command.replace('\n', ' '))

    p.write_output(extract_task.out, output_root)
    if gene_map_ht_path:
        p.write_output(extract_task.group_file, f'{output_root}.gene.txt')
    p.write_output(extract_task.stdout, f'{output_root}.log')
    return extract_task


def export_pheno(p: Pipeline, output_path: str, pheno: str, coding: str, trait_type: str, module: str,
                 docker_image: str, n_threads: int = 8, storage: str = '500Mi', additional_args: str = ''):
    extract_task: pipeline.pipeline.Task = p.new_task(name='extract_pheno',
                                                      attributes={
                                                          'pheno': pheno,
                                                          'coding': coding,
                                                          'trait_type': trait_type
                                                      })
    extract_task.image(docker_image).cpu(n_threads).storage(storage)
    python_command = f"""set -o pipefail; python3 {SCRIPT_DIR}/export_pheno.py
    --load_module {module}
    --trait_type {trait_type}
    --pheno {shq(pheno)} --sex both_sexes
    {"--additional_args " + additional_args if additional_args else ''}
    {"--coding " + coding if coding else ''}
    --output_file {extract_task.out}
    --n_threads {n_threads} | tee {extract_task.stdout}
    ; """.replace('\n', ' ')

    extract_task.command(python_command)

    p.write_output(extract_task.out, output_path)
    p.write_output(extract_task.stdout, f'{output_path}.log')
    return extract_task


def fit_null_glmm(p: Pipeline, output_root: str, pheno_file: pipeline.pipeline.Resource, trait_type: str, covariates: str,
                  plink_file_root: str, docker_image: str, sparse_grm: pipeline.pipeline.Resource = None,
                  sparse_grm_extension: str = None, skip_model_fitting: bool = False,
                  n_threads: int = 16, storage: str = '1500Mi', memory: str = '60G'):
    analysis_type = "variant" if sparse_grm is None else "gene"
    pheno_col = 'value'
    user_id_col = 'userId'
    in_bfile = p.read_input_group(**{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    fit_null_task = p.new_task(name=f'fit_null_model',
                               attributes={
                                   'analysis_type': analysis_type,
                                   'trait_type': trait_type
                               }).cpu(n_threads).storage(storage).image(docker_image).memory(memory)
    output_files = {ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}' for ext in
                   ('rda', '_30markers.SAIGE.results.txt', f'{analysis_type}.varianceRatio.txt')}
    if analysis_type == 'gene':
        sparse_sigma_extension = sparse_grm_extension.replace("GRM", "Sigma")
        output_files[f'{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'] = \
            f'{{root}}.{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'
    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f'perl -pi -e s/^chr// {in_bfile.bim}'
    # if trait_type == 'icd':
    #     bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
    #                         f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")

    command = (f'set -o pipefail; Rscript /usr/local/bin/step1_fitNULLGLMM.R '
               f'--plinkFile={in_bfile} '
               f'--phenoFile={pheno_file} '
               f'--covarColList={covariates} '
               f'--phenoCol={pheno_col} '
               f'--sampleIDColinphenoFile={user_id_col} '
               f'--traitType={saige_pheno_types[trait_type]} '
               f'--outputPrefix={fit_null_task.null_glmm} '
               f'--outputPrefix_varRatio={fit_null_task.null_glmm}.{analysis_type} '
               f'--skipModelFitting={str(skip_model_fitting).upper()} ')
    if analysis_type == "gene":
        fit_null_task.declare_resource_group(sparse_sigma={sparse_sigma_extension: f'{{root}}.{sparse_sigma_extension}'})
        command += (f'--IsSparseKin=TRUE '
                    f'--sparseGRMFile={sparse_grm[sparse_grm_extension]} '
                    f'--sparseGRMSampleIDFile={sparse_grm[f"{sparse_grm_extension}.sampleIDs.txt"]} '
                    f'--isCateVarianceRatio=TRUE ')
    command += f'--nThreads={n_threads} --LOCO=FALSE 2>&1 | tee {fit_null_task.stdout}'
    command = '; '.join([bim_fix_command, command])
    fit_null_task.command(command)
    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f'{output_root}.{analysis_type}.log')
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task


def run_saige(p: Pipeline, output_root: str, model_file: str, variance_ratio_file: str,
              vcf_file: pipeline.pipeline.ResourceGroup, samples_file: pipeline.pipeline.ResourceGroup,
              docker_image: str,
              group_file: str = None, sparse_sigma_file: str = None, use_bgen: bool = True,
              trait_type: str = 'continuous',
              chrom: str = 'chr1', min_mac: int = 1, min_maf: float = 0, max_maf: float = 0.5,
              memory: str = '', storage: str = '500Mi'):

    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    run_saige_task: pipeline.pipeline.Task = p.new_task(name=f'run_saige',
                                                        attributes={
                                                            'analysis_type': analysis_type
                                                        }).cpu(1).storage(storage).image(docker_image)  # Step 2 is single-threaded only

    if analysis_type == 'gene':
        run_saige_task.declare_resource_group(result={'gene.txt': '{root}',
                                                      'single.txt': '{root}_single'})
    else:
        run_saige_task.declare_resource_group(result={'single_variant.txt': '{root}'})

    command = (f'set -o pipefail; {MKL_OFF} Rscript /usr/local/bin/step2_SPAtests.R '
               f'--minMAF={min_maf} '
               f'--minMAC={min_mac} '
               f'--maxMAFforGroupTest={max_maf} '
               f'--sampleFile={samples_file} '
               f'--GMMATmodelFile={model_file} '
               f'--varianceRatioFile={variance_ratio_file} '
               f'--SAIGEOutputFile={run_saige_task.result} ')

    if use_bgen:
        command += (f'--bgenFile={vcf_file.bgen} '
                    f'--bgenFileIndex={vcf_file["bgen.bgi"]} ')
    else:
        command += (f'--vcfFile={vcf_file["vcf.gz"]} '
                    f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
                    f'--chrom={chrom} '
                    f'--vcfField=GT ')
    if analysis_type == "gene":
        if saige_pheno_types[trait_type] == 'binary':
            command += f'--IsOutputPvalueNAinGroupTestforBinary=TRUE '
        command += (f'--groupFile={group_file} '
                    f'--sparseSigmaFile={sparse_sigma_file} '
                    f'--IsSingleVarinGroupTest=TRUE '
                    f'--IsOutputBETASEinBurdenTest=TRUE ')
    command += f'--IsOutputAFinCaseCtrl=TRUE 2>&1 | tee {run_saige_task.stdout}; '
    if analysis_type == 'gene':
        command += f"input_length=$(wc -l {group_file} | awk '{{print $1}}'); " \
            f"output_length=$(wc -l {run_saige_task.result['gene.txt']} | awk '{{print $1}}'); " \
            f"echo 'Got input:' $input_length 'output:' $output_length | tee -a {run_saige_task.stdout}; " \
            f"if [[ $input_length > 0 ]]; then echo 'got input' | tee -a {run_saige_task.stdout}; " \
            f"if [[ $output_length == 1 ]]; then echo 'but not enough output' | tee -a {run_saige_task.stdout}; exit 1; fi; fi"
    run_saige_task.command(command)
    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f'{output_root}.{analysis_type}.log')
    return run_saige_task


def load_results_into_hail(p: Pipeline, output_root: str, pheno: str, coding: str, trait_type: str, tasks_to_hold,
                           vep_path: str, docker_image: str, gene_map_path: str = None, reference: str = 'GRCh38',
                           analysis_type: str = 'gene', n_threads: int = 8, storage: str = '500Mi'):

    load_data_task: pipeline.pipeline.Task = p.new_task(name=f'load_data',
                                                        attributes={
                                                            'pheno': pheno,
                                                            'coding': coding,
                                                            'trait_type': trait_type
                                                        }).image(docker_image).cpu(n_threads).storage(storage)
    load_data_task.always_run().depends_on(*tasks_to_hold)
    python_command = f"""python3 {SCRIPT_DIR}/load_results.py
    --input_dir {shq(output_root)}
    --pheno {shq(pheno)}
    {"--coding " + coding if coding else ''}
    --trait_type {trait_type}
    {"--gene_map_ht_raw_path " + gene_map_path if gene_map_path else ''}
    --ukb_vep_ht_path {vep_path}
    --overwrite --reference {reference}
    --analysis_type {analysis_type}
    --n_threads {n_threads} | tee {load_data_task.stdout}
    ;""".replace('\n', ' ')

    python_command = python_command.replace('\n', '; ').strip()
    command = f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" ' + python_command
    load_data_task.command(command)
    p.write_output(load_data_task.stdout, f'{output_root}/{pheno}_loading.log')
    return load_data_task


def qq_plot_results(p: Pipeline, output_root: str, tasks_to_hold, export_docker_image: str, R_docker_image: str,
                    n_threads: int = 8, storage: str = '500Mi'):

    qq_export_task: pipeline.pipeline.Task = p.new_task(name='qq_export').image(export_docker_image).cpu(n_threads).storage(storage)
    qq_export_task.always_run().depends_on(*tasks_to_hold)

    python_command = f"""python3 {SCRIPT_DIR}/export_results_for_qq.py
    --input_dir {shq(output_root)}
    --output_file {qq_export_task.out}
    --n_threads {n_threads}
    ; """.replace('\n', ' ')

    command = f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" ' + python_command
    qq_export_task.command(command)

    qq_task: pipeline.pipeline.Task = p.new_task(name='qq_plot').image(R_docker_image).cpu(n_threads).storage(storage).always_run()
    qq_task.declare_resource_group(result={ext: f'{{root}}_Pvalue_{ext}'
                                           for ext in ('qqplot.png', 'manhattan.png', 'manhattan_loglog.png', 'qquantiles.txt')})
    R_command = f"/saige-pipelines/scripts/qqplot.R -f {qq_export_task.out} -o {qq_task.result} -p Pvalue; "
    qq_task.command(R_command)

    p.write_output(qq_task.result, output_root)
    return qq_export_task, qq_task


def get_tasks_from_pipeline(p):
    return dict(Counter(map(lambda x: x.name, p.select_tasks(""))))


def load_jobs_by_batch_ids(batch_ids, billing_project: str = 'ukb_diverse_pops'):
    import batch_client
    bc = batch_client.client.BatchClient(billing_project=billing_project)
    if isinstance(batch_ids, int):
        batch_ids = [batch_ids]
    all_jobs = []
    for batch_id in batch_ids:
        batch = bc.get_batch(batch_id)
        all_jobs.extend(list(batch.jobs()))
    return all_jobs


def get_costs_by_attribute(attributes, jobs = None, batch_ids = None, filter_job_name: str = None, get_status_instead: bool = False,
                           billing_project: str = 'ukb_diverse_pops'):
    import warnings
    if batch_ids is not None:
        jobs = load_jobs_by_batch_ids(batch_ids, billing_project=billing_project)
    if jobs is None:
        warnings.warn('batch_ids and jobs are both None. Bailing...')
        return None

    from collections import defaultdict
    summary = defaultdict(float) if not get_status_instead else defaultdict(lambda: defaultdict(int))
    for job in jobs:
        if filter_job_name is not None and job['attributes'].get('name') != filter_job_name:
            continue
        key = tuple(job['attributes'].get(attribute, '') for attribute in attributes)
        if get_status_instead:
            summary[key][job['state']] += 1
        else:
            summary[key] += float(job['cost'].lstrip('$'))

    ret_summary = {}
    for k, v in summary.items():
        if get_status_instead:
            if k not in ret_summary: ret_summary[k] = {}
            for k2, v2 in v.items():
                ret_summary[k][k2] = v2
        else:
            ret_summary[k] = round(v, 4)
    return ret_summary
