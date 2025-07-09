
include { TIgGERBayesianGenotypeInference } from './tigger/main.nf'
include { PIgLETGenotypeInference } from './piglet/main.nf'

workflow GenotypeInference {

  take:
    genotype_input
    segment
    method
    find_unmutated
    single_assignments
    threshold_file
    z_threshold   

  main:
    if (method == "piglet") {
      PIgLETGenotypeInference(
        genotype_input,
        threshold_file,
        segment,
        find_unmutated,
        single_assignments,
        z_threshold
      ).set { results }
    } else if (method == "tigger_bayesian") {
      TIgGERBayesianGenotypeInference(
        genotype_input,
        segment,
        find_unmutated,
        single_assignments
      ).set { results }
    } else {
      error "Unsupported genotype method: ${method}"
    }

  emit:
    genotype = results.genotype
    reference = results.reference
}
