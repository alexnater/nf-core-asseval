//
// This file holds several functions specific to the workflow/asseval.nf in the nf-core/asseval pipeline
//

import groovy.json.JsonSlurper
import groovy.yaml.YamlBuilder

class WorkflowAssEval {

    //
    // Function to parse FastP json files to get total number of reads before and after trimming per sample
    //
    public static Tuple getNumberOfReads(json_files) {
        def total_reads = 0
        def trimmed_reads = 0
        json_files.each { json_file ->
            def jsonSlurper = new JsonSlurper()
            def json = jsonSlurper.parse(json_file)
            total_reads += json.summary.before_filtering.total_reads
            trimmed_reads += json.summary.after_filtering.total_reads
        }
        return new Tuple(total_reads, trimmed_reads)
    }

    //
    // Function to parse fasta index files and get a list of contigs in the assembly
    //
    public static List getContigs(fai_file, min_size=0) {
        def contigs = []
        fai_file.eachLine { line ->
            def fields = line.trim().split('\t')
            if (fields[1] as int >= min_size) {
                contigs << fields[0]
            }
        }
        return contigs
    }

    //
    // Function to parse fasta index files and get a list of regions corresponding to all contigs in the assembly
    //
    public static List getContigRegions(fai_file, min_size=0) {
        def regions = []
        fai_file.eachLine { line ->
            def fields = line.trim().split('\t')
            if (fields[1] as int >= min_size) {
                regions << "${fields[0]}:1-${fields[1]}"
            }
        }
        return regions
    }

    //
    // Function to generate YAML file
    //
    public static void printYAML(outfile, hifi_depth=0, ont_depth=0, hic_depth=0) {
        def builder = new YamlBuilder()
        builder.depth {
            hifi "${hifi_depth}"
            ont "${ont_depth}"
            hic "${hic_depth}"
        }
        println(builder.toString())
    }
}