#!/usr/bin/env bash

#------------------------------------------------------------------------------#

## REQUIRED: path to clinvar_workflow package
PKG_PATH=/Users/elisabethmlynarski/repos/bitbucket/cv_wkflw_dev/

################## FileIO ######################################################
## REQUIRED: input variant FILE - relative or absolute path
VAR_FILE=demo_input_variant_files/demo_variants_ASD_hg19.txt

## REQUIRED: output file DIRECTORY - relative or absolute path
OUT_DIR=demo_output

## REQUIRED: prefix for the outputs (default = '')
OUT_PREFIX='demo_script_ASD'

################## update values based on YOUR current input file ##############
## REQUIRED: Genome build: hg19 or hg38 (default = hg19)
BUILD='hg19'

## REQUIRED: 4 variant columns (names can change but order must remain the same!)
COLS_VAR='CHR,POS,REF,ALT'

## optional: list of additional input columns to include in output DF
COLS_INPUT='VAR_FUNCTION,VAR_CATEGORY,STUDY_TYPE,VAR_VALIDATION'

#------------------------------------------------------------------------------#

PY_SCRIPT="${PKG_PATH}"clinvar_workflow_scripts/clinvar_exploratory_analysis.py


python $PY_SCRIPT --pkg_path $PKG_PATH --var_file $VAR_FILE --out_dir $OUT_DIR --out_prefix $OUT_PREFIX --build $BUILD --cols_var $COLS_VAR --cols_input "${COLS_INPUT}"


