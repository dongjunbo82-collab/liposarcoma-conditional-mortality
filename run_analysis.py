"""Run locally from a code checkout; restricted inputs are never uploaded."""
from pathlib import Path
import argparse,os,subprocess,sys
ROOT=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['primary','secondary','models','smr','figures','all'],default='primary')
    args=parser.parse_args();os.chdir(ROOT)
    base=Path(os.environ.get('LPS_OUTPUT_DIR','outputs/population_study'))
    for suffix in ['','revision','survivorship','revision/focus','revision/cancer','revision/smr','standardized']:(base/suffix).mkdir(parents=True,exist_ok=True)
    r=os.environ.get('RSCRIPT','Rscript')
    def R(name,*extra):subprocess.run([r,'analysis/'+name,*extra],check=True)
    def py(name):subprocess.run([sys.executable,'analysis/'+name],check=True)
    R('check_dependencies.R')
    if args.stage in ['primary','all']:
        raw=Path(os.environ.get('SEER_EXPORT','inputs/seer_export.tsv'))
        if not raw.is_file():raise SystemExit('Provide your authorized local 23-column SEER export; see docs/input_schema.md.')
        R('population_primary.R');py('verify_primary_estimates.py');py('followup_support.py')
    if args.stage in ['secondary','all']:
        py('survivorship_cod_map.py');R('survivorship_extension.R');R('revision_analyses.R');R('clinical_context_extension.R');R('cancer_period_sensitivity.R')
    if args.stage in ['models','all']:
        R('population_standardized.R','500');R('age_definition_sensitivity.R');R('age_model_checks.R')
    if args.stage in ['smr','all']:
        py('cardiovascular_smr.py');R('validate_cardiovascular_smr.R')
    if args.stage in ['figures','all']:R('focus_figures.R')

if __name__=='__main__':main()
