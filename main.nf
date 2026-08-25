#!/usr/bin/env nextflow
// RiboFlow (genome) — DSL2 entrypoint.
nextflow.enable.dsl = 2

include { RIBOFLOW } from './workflows/riboflow.nf'

workflow {
    RIBOFLOW()
}
