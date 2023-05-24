version 1.0

task RunVarpipeline {
  input {
    File R1
    File R2
    File REF
    String SAMPL
  }
  command {
    Varpipeline -q ${R1} -q2 ${R2} -r ${REF} -g NC_000962 -n ${SAMPL} -v -a -k
  }
  runtime {
    docker: "dbest/varpipe3:latest"
  }
}

workflow WfVarPipeline {
  input {
    File R1
    File R2
    File REF
    String SAMPL
  }
  call RunVarpipeline {
    input:
    R1 = R1,
    R2 = R2,
    REF = REF,
    SAMPL = SAMPL
  }
}
