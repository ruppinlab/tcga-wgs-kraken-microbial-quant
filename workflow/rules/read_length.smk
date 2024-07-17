rule bbmap_read_length_histogram_pe:
    input:
        fq1=BOWTIE2_FILTERED_FASTQ_R1_FILE,
        fq2=BOWTIE2_FILTERED_FASTQ_R2_FILE,
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


rule bbmap_read_length_histogram_se:
    input:
        fq1=BOWTIE2_FILTERED_FASTQ_SE_FILE,
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


rule bbmap_max_read_length:
    input:
        READ_LENGTH_HISTOGRAM_FILE,
    output:
        READ_LENGTH_FILE,
    log:
        READ_LENGTH_LOG,
    run:
        import pandas as pd

        cols = [
            "length",
            "reads",
            "pct_reads",
            "cum_reads",
            "cum_pct_reads",
            "bases",
            "pct_bases",
            "cum_bases",
            "cum_pct_bases",
        ]
        print(f"Getting maximum read length from {input[0]}", flush=True)
        df = pd.read_csv(input[0], sep="\t", comment="#", names=cols)
        max_length = df["length"].max().astype(str)
        with open(output[0], "w") as fh:
            fh.write(max_length)
        with open(log[0], "w") as fh:
            fh.write(f"Wrote max length {max_length} file")
