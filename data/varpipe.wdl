task RunVarpipeline {
  command {
    docker run --rm -v ${DATA}:${DAT} -w ${TOOLS} dbest/varpipe:latest Varpipeline -q ${R1} -q2 ${R2} -r ${REF} -g NC_000962 -n ${SAMPL} -v -a -k
  }
  
  runtime {
    docker: "dbest/varpipe3:latest"
  }
  
  input {
    File DATA
    String DAT
    String TOOLS
    File R1
    File R2
    File REF
    String SAMPL
  }
}

workflow MyVarPipeline {
  call RunVarpipeline {
    input:
      DATA = "path/to/data",
      DAT = "$DAT",
      TOOLS = "$TOOLS",
      R1 = file("path/to/R1.fastq"),
      R2 = file("path/to/R2.fastq"),
      REF = file("path/to/reference.fa"),
      SAMPL = "$SAMPL"
  }
}
