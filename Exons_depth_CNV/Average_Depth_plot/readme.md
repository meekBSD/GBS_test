
    java -jar ${picard} CollectHsMetrics \
        R=${hg19} \
        I=${sample}.consensus.filter.BAM \
        O=${sample}.target.stat.xls \
        BI=${intervalF} \
        TI=${intervalF} \
        PER_TARGET_COVERAGE=${sample}.target.coverage \
        PER_BASE_COVERAGE=${sample}.base.coverage \
        MQ=0 \
        Q=0 \
        COVMAX=100000 \
        CLIP_OVERLAPPING_READS=true

input: bed file and base coverage

