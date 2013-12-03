package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE variant calling pipeline with just the UnifiedGenotyper, VariantEval and optional VariantFiltration.
 * For a complete description of the suggested for a variant calling pipeline see the latest version of the Best Practice Variant Detection document
 */
class ExampleUnifiedGenotyper extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Input(doc="Rod file.", shortName="D", required=false)
  var dbsnpFile: File = _

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    // Set the memory limit to 8 gigabytes on each command.
    this.memoryLimit = 8
  }

  def script() {
    // Create the four functions that we may run depending on options.
    val genotyper = new UnifiedGenotyper with UnifiedGenotyperArguments

    genotyper.scatterCount = 10
    genotyper.input_file :+= qscript.bamFile
    genotyper.out = swapExt(qscript.bamFile, "clean.dedup.recal.bam", "raw.vcf")
    genotyper.dbsnp = qscript.dbsnpFile
    genotyper.logging_level = "INFO"
    genotyper.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    genotyper.output_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY

    add(genotyper)
  }
}
