#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}
if (params.workflow == "rnafusion"){
  include {RNAFUSION} from './workflows/rnafusion'
}
if (params.workflow == "cnv"){
  include {CNV} from './workflows/cnv'
}
// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "rnaseq"){
    RNASEQ()
    }
  if (params.workflow == "rnafusion"){
    RNAFUSION()
    }
  if (params.workflow == "cnv"){
    CNV()
  }
}
