#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}
// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "rnaseq"){
    RNASEQ()
    }
}
