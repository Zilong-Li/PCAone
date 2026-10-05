#!/usr/bin/env nextflow

/*
 * LD decay curves, standard vs ancestry-adjusted, for every population in
 * params.data. See docs/guide/plotting.md ("LD decay").
 *
 *   make -C scripts summarise_ld_r2bin      # once, builds the binning tool
 *   nextflow run workflows/ld.nf --data data --pops giraffe --K 5,10 --run_step plot
 *
 * run_step
 *   curve  LD decay curves             -> ld_curve_k<K>.png
 *   cross  mean r2 of unlinked pairs   -> cross_mean_{adj_<K>,std}.txt
 *   plot   both, in one plot           -> ld_combine_k<K>.png
 *
 * Needs PCAone and plink (1.9) in PATH, and Rscript with data.table.
 */

params.run_step = "curve"
params.pops = ["giraffe"]
params.data = "data"                          // <pop>.{bed,bim,fam}
params.scripts = "${projectDir}/../scripts"
params.results = "results"
params.maf = 0.05
params.thin = 0.05
params.ld_bp = 1000000
params.ld_bins = 40                           // log-spaced distance bins
params.K = [10]
params.perm_chr = "1,2"                       // chromosomes shuffled for the unlinked pairs
params.correct = true                         // correct r2 for sample size

log.info """\
         L D   D E C A Y
         ===================================
      run_step:       : ${params.run_step}
          pops:       : ${params.pops}
           maf:       : ${params.maf}
          thin:       : ${params.thin}
        window:       : ${params.ld_bp}
             K:       : ${params.K}
         """
         .stripIndent()

process qc_filters {
    input:
    tuple val(pop), path(bfiles)

    output:
    tuple val(pop),
        path("${pop}.bed"),
        path("${pop}.bim"),
        path("${pop}.fam")

    script:
    base = bfiles[0].baseName
    """
    plink --bfile $base --maf ${params.maf} --max-maf 0.499 --thin ${params.thin} --make-bed --out $pop
    """
}

// Unlinked pairs: give the SNPs of params.perm_chr random chromosome labels.
// plink re-sorts them, so an LD window now pairs SNPs from different original
// chromosomes. The SNP id keeps the original chromosome ("<chr>:<bp>").
process permute_plink {
    input:
    tuple val(pop), path(bed), path(bim), path(fam)

    output:
    tuple val(pop),
        path("${pop}.perm.bed"),
        path("${pop}.perm.bim"),
        path("${pop}.perm.fam")

    script:
    base = bed.baseName
    """
    plink --bfile $base --chr ${params.perm_chr} --make-bed --out ${pop}.chr
    R -s -e 'd=read.table("${pop}.chr.bim",colClasses="character"); d[,2]=paste0(d[,1],":",d[,4]); d[,1]=sample(d[,1]); write.table(d,file="${pop}.chr.bim.bk",sep="\\t",quote=F,col.names=F,row.names=F)'
    mv ${pop}.chr.bim.bk ${pop}.chr.bim
    plink --bfile ${pop}.chr --make-bed --out ${pop}.perm
    """
}

process adj_ld_matrix {
    input:
    tuple val(K), val(pop), path(bed), path(bim), path(fam)

    output:
    tuple val(pop), val(K), path(bed), path(bim), path(fam), path("adj.${K}.eigvecs")

    script:
    base = bed.baseName
    // --scale 0 gives the unstandardized PCs the removed -D/--ld used
    """
    PCAone --bfile $base -k ${K} -d 0 --scale 0 --out adj.${K}
    """
}

process adj_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}"

    input:
    tuple val(pop), val(K), path(bed), path(bim), path(fam), path(eigvecs)

    output:
    tuple val(pop), val(K), path("adj.${K}.ld.gz")

    script:
    base = bed.baseName
    """
    PCAone --bfile $base --USV adj.${K} --ld-bp ${params.ld_bp} --print-r2 --out adj.${K}
    """
}

process std_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", pattern: "*.ld.gz"

    input:
    tuple val(pop), path(bed), path(bim), path(fam)

    output:
    tuple val(pop), path("std.ld.gz"), path(fam)

    script:
    base = bed.baseName
    """
    PCAone --bfile $base --ld-stats 1 --ld-bp ${params.ld_bp} --print-r2 --out std
    """
}

// mean r2 of the pairs whose SNPs come from different original chromosomes,
// with the PCs of the whole data set for the adjusted LD
process cross_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'

    input:
    tuple val(pop), val(K), path(bed), path(bim), path(fam), path(eigvecs)

    output:
    tuple val(pop), val(K), path("cross_mean_adj_${K}.txt"), path("cross_mean_std.txt")

    script:
    base = bed.baseName
    """
    PCAone --bfile $base --USV adj.${K} --ld-bp ${params.ld_bp} --print-r2 --out perm_adj
    PCAone --bfile $base --ld-stats 1 --ld-bp ${params.ld_bp} --print-r2 --out perm_std
    for t in adj std; do
      zcat perm_\$t.ld.gz | awk 'NR > 1 { split(\$3, a, ":"); split(\$6, b, ":"); if (a[1] != b[1]) { s += \$7; n++ } }
        END { if (n == 0) { print "no unlinked pairs" > "/dev/stderr"; exit 1 }; printf "%.7g\\n", s / n }' > cross_mean_\$t.txt
    done
    mv cross_mean_adj.txt cross_mean_adj_${K}.txt
    """
}

process make_adj_ld_bin {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true

    input:
    tuple val(pop), val(K), path(adjr2)

    output:
    tuple val(pop), val(K), path("adj.${K}.decay.tsv")

    script:
    """
    ${params.scripts}/summarise_ld_r2bin -i ${adjr2} --max ${params.ld_bp} --bins ${params.ld_bins} -o adj.${K}.decay.tsv
    """
}

process make_std_ld_bin {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, pattern: "*.decay.tsv"

    input:
    tuple val(pop), path(stdr2), path(fam)

    output:
    tuple val(pop), path("std.decay.tsv"), path(fam)

    script:
    """
    ${params.scripts}/summarise_ld_r2bin -i ${stdr2} --max ${params.ld_bp} --bins ${params.ld_bins} -o std.decay.tsv
    """
}

process plot_ld_curve {
    cache false
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'

    input:
    tuple val(pop), val(K), path(adj), path(std), path(fam)

    output:
    path("ld_curve_k${K}.png")

    script:
    correct = params.correct ? "--correct" : ""
    """
    Rscript ${params.scripts}/plot-ld-decay.R ${adj} ${std} --labels Adjusted,Standard \
        -n ${fam} ${correct} --title "${pop}, K=${K}" -o ld_curve_k${K}.png
    """
}

process plot_ld_combine {
    cache false
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'

    input:
    tuple val(pop), val(K), path(adj), path(std), path(fam), path(cross_adj), path(cross_std)

    output:
    path("ld_combine_k${K}.png")

    script:
    correct = params.correct ? "--correct" : ""
    """
    Rscript ${params.scripts}/plot-ld-decay.R ${adj} ${std} --labels Adjusted,Standard \
        -n ${fam} ${correct} --baseline \$(cat ${cross_adj}),\$(cat ${cross_std}) \
        --title "${pop}, K=${K}" -o ld_combine_k${K}.png
    """
}

// --pops a,b and --K 5,10 on the command line, or lists in a config
def as_list(x) { x instanceof List ? x : x.toString().tokenize(',')*.trim() }

// pcs: pop, K, bed, bim, fam, eigvecs
workflow ld_curve {
    take:
    data
    pcs

    main:
    adj = adj_ld_r2(pcs) | make_adj_ld_bin
    std = std_ld_r2(data) | make_std_ld_bin
    ld = adj.combine(std, by: 0)               // pop, K, adj, std, fam
    plot_ld_curve(ld)

    emit:
    ld
}

workflow ld_cross {
    take:
    data
    pcs

    main:
    perm = permute_plink(data)
    ld = perm.combine(pcs.map { pop, K, bed, bim, fam, eig -> tuple(pop, K, eig) }, by: 0)
        .map { pop, bed, bim, fam, K, eig -> tuple(pop, K, bed, bim, fam, eig) } | cross_ld_r2

    emit:
    ld                                        // pop, K, cross_adj, cross_std
}

workflow {
    ch_plink = channel.fromFilePairs("${params.data}/*.{bed,bim,fam}", size:3, checkIfExists:true) {
        file -> file.baseName
    }.filter {key, files -> key in as_list(params.pops)}

    ch_K = channel.fromList(as_list(params.K))
    ch_data = qc_filters(ch_plink)
    ch_pcs = adj_ld_matrix(ch_K.combine(ch_data))

    if( params.run_step == 'curve') {
        ld_curve(ch_data, ch_pcs)
    }
    if( params.run_step == 'cross') {
        ld_cross(ch_data, ch_pcs)
    }
    if( params.run_step == 'plot') {
        curve = ld_curve(ch_data, ch_pcs)
        cross = ld_cross(ch_data, ch_pcs)
        combine = curve.combine(cross, by: [0,1])
        plot_ld_combine(combine)
    }
}
