#!/usr/bin/env nextflow


/*
 * Define the default parameters
 */

params.run_step = "curve"
params.pops = ["giraffe"]
params.data = "data"
params.scripts = "scripts"
params.results = "results"
params.maf=0.05
params.thin=0.05
params.ld_bp=1000000
params.K=[10]
params.py27="/maps/projects/alab/people/rlk420/software/miniforge3/envs/py27/bin/python2"

log.info """\
         R N A S E Q - N F   P I P E L I N E
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
    plink --bfile $base --chr 1,2 --make-bed --out ${pop}.chr
    R -s -e 'd=read.table("${pop}.chr.bim"); d[,1]=sample(d[,1]);write.table(d,file="${pop}.chr.bim.bk",sep="\\t",quote=F,col.names=F,row.names=F)'
    mv ${pop}.chr.bim.bk ${pop}.chr.bim
    plink --bfile ${pop}.chr --make-bed --out ${pop}.perm
    """
}

process adj_ld_matrix {
    // publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}"
    
    input:
    tuple val(K), val(pop), path(bed), path(bim), path(fam)

    output:
    tuple val(pop), val(K), path("adj.${K}.residuals"), emit: residuals
    tuple val(pop), val(K), path("adj.${K}.kept.bim"), emit: kept

    script:
    base = bed.baseName
    """
    PCAone --bfile $base -k ${K} -d 0 --ld-stats 0 --ld --out adj.${K}
    """
}

process adj_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}"
    
    input:
    tuple val(pop), val(K), path(residuals)
    tuple val(pop), val(K), path(kept)

    output:
    tuple val(pop), val(K), path("adj.${K}.ld.gz"), emit: adjr2

    script:
    """
    PCAone -B ${residuals} --ld-bim ${kept} --ld-bp ${params.ld_bp} --print-r2 --out adj.${K}
    """
}

process std_ld_matrix {
    // publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}"
    
    input:
    tuple val(pop), path(bed), path(bim), path(fam)

    output:
    tuple val(pop), path("std.residuals"), path("std.kept.bim")

    script:
    base = bed.baseName
    """
    PCAone --bfile $base -k 1 -d 0 --ld-stats 1 --ld --out std 
    """
}

process std_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}"
    
    input:
    tuple val(pop), path("std.residuals"), path("std.kept.bim")

    output:
    tuple val(pop), path("std.ld.gz"), emit: stdr2

    script:
    """
    PCAone -B std.residuals --ld-bim std.kept.bim --ld-bp ${params.ld_bp} --print-r2 --out std
    """
}

process cross_std_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'
    
    input:
    tuple val(pop), path(stdr2)

    output:
    tuple val(pop), path("cross_mean_std.txt"), emit: cross_std
    
    script:
    """
    zcat ${stdr2} | grep '1:' | grep '2:' | awk '{s+=\$7}END{print s/NR}' > cross_mean_std.txt
    """
}

process cross_adj_ld_r2 {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'
    
    input:
    tuple val(pop), val(K), path(adjr2)

    output:
    tuple val(pop), val(K), path("cross_mean_adj_${K}.txt"), emit:cross_adj
    
    script:
    """
    zcat ${adjr2} | grep '1:' | grep '2:' | awk '{s+=\$7}END{print s/NR}' > cross_mean_adj_${K}.txt
    """
}

process make_adj_ld_bin {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true

    input:
    tuple val(pop), val(K), path(adjr2)

    output:
    tuple val(pop), val(K), path("adj.${K}.decay_bins"), emit: adj

    script:
    """
    ${params.py27} /maps/projects/alab/people/rlk420/ld/paper/ld_decay_calc.py -i ${adjr2} -o adj.${K}
    """
}

process make_std_ld_bin {
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true
    
    input:
    tuple val(pop), path(stdr2)

    output:
    tuple val(pop), path("std.decay_bins"), emit: std

    script:
    """
    ${params.py27} /maps/projects/alab/people/rlk420/ld/paper/ld_decay_calc.py -i ${stdr2} -o std
    """
}

process plot_ld_curve {
    cache false
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'

    input:
    tuple val(pop), val(K), path(adj), path(std)
    
    output:
    path("ld_curve_k${K}.png")
    
    script:
    """
    #!/usr/bin/env Rscript
    source("${params.scripts}/common.R")
    bitmap("ld_curve_k${K}.png", res=300, w=4,h=4)
    pops <- c("giraffe"=29, "ASW-CEU-YRI"=150) ## number of samples
    adj <- parse_ld_bins("${adj}", n = pops["${pop}"])
    std <- parse_ld_bins("${std}", n = pops["${pop}"])
    lwd <- 4
    plot_ld_curve(adj, lwd=lwd, main = paste0("${pop}, K=${K}, mean=", round(mean(adj[,"r2"], na.rm=T), 3)), cex.main=2)
    lines(std[,"dist"], std[,"r2"], lty = 2, lwd = lwd, col=2)
    legend("topright",lty=c(2,1),c("Standard","Adjusted"),lwd=lwd,bty="n", cex = 2)
    dev.off()
    """
}

process plot_ld_combine {
    cache false
    publishDir "${params.results}/${pop}/thin_${params.thin}/${params.run_step}", overwrite:true, mode:'copy'

    input:
    tuple val(pop), val(K), path(adj), path(std), path(cross_adj), path(cross_std)

    output:
    path("ld_combine_k${K}.png")
    
    script:
    """
    #!/usr/bin/env Rscript
    source("${params.scripts}/common.R")
    bitmap("ld_combine_k${K}.png", res=300, w=4,h=4)
    pops <- c("giraffe"=29, "ASW-CEU-YRI"=150) ## number of samples
    adj <- parse_ld_bins("${adj}", n = pops["${pop}"])
    std <- parse_ld_bins("${std}", n = pops["${pop}"])
    a <- scan("${cross_adj}")
    b <- scan("${cross_std}")
    cross <- c(bs(a, n = pops["${pop}"]), bs(b, n = pops["${pop}"]))
    lwd <- 4
    par(mar = c(6, 5, 5, 6))
    plot_ld_curve(adj, lwd=lwd, main = paste0("${pop}, K=${K}, mean=", round(mean(adj[,"r2"], na.rm=T), 3)), cex.main=2)
    lines(std, lty = 2, lwd = lwd, col=2)
    a <- cross[1]
    b <- cross[2]
    dist <- adj[,1]
    lines(c(tail(dist,1)*c(1.4,2.5)),rep(a, 2),lty=1,lwd=lwd,xpd=T, col=1)
    lines(c(tail(dist,1)*c(1.4,2.5)),rep(b, 2),lty=2,lwd=lwd,xpd=T, col=2)
    mtext("Diff\nChr",las=1,at=2*5e6,xpd=T,side=1,font=1,line=2,cex = 1.2)
    legend("topright",c("Standard","Adjusted"),lty=c(2,1), col = c(2,1),lwd=lwd,bty="n", cex = 2)
    dev.off()
    """
}

workflow ld_curve {
    take:
    data
    K

    main:
    adj = adj_ld_matrix(K.combine(data)) | adj_ld_r2 | make_adj_ld_bin
    std = std_ld_matrix(data) | std_ld_r2 | make_std_ld_bin
    ld = adj.combine(std, by: 0)
    // ld.view()
    plot_ld_curve(ld)
    
    emit:
    ld
}

workflow ld_cross {
    take:
    data
    K

    main:
    ch_perm = permute_plink(data)
    adj = adj_ld_matrix(K.combine(ch_perm)) | adj_ld_r2 | cross_adj_ld_r2
    std = std_ld_matrix(ch_perm) | std_ld_r2 | cross_std_ld_r2
    ld = adj.combine(std, by: 0)

    emit:
    ld
}

workflow {
    ch_plink = channel.fromFilePairs("${params.data}/*.{bed,bim,fam}", size:3, checkIfExists:true) {
        file -> file.baseName
    }.filter {key, files -> key in params.pops}

    ch_K = channel.fromList(params.K)
    ch_data = qc_filters(ch_plink)

    if( params.run_step == 'curve') {
        ld_curve(ch_data, ch_K)
    }
    if( params.run_step == 'cross') {
        ld_cross(ch_data, ch_K)
    }
    if( params.run_step == 'plot') {
        curve = ld_curve(ch_data, ch_K)
        cross = ld_cross(ch_data, ch_K)
        combine = curve.combine(cross, by: [0,1])
        combine.view()
        plot_ld_combine(combine)
    }
}
