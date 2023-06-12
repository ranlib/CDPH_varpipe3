version 1.0

task RunVarpipeline {
  input {
    File R1
    File R2
    File REF
    String SAMPL
    Int threads = 1
    String genome = "NC_000962"
    String memory = "10GB"
  }
  command {
    Varpipeline -q ${R1} -q2 ${R2} -r ${REF} -g ${genome} -n ${SAMPL} -t ${threads} -v -a -k
  }
  runtime {
    docker: "dbest/varpipe3:latest"
    cpu: threads
    memory: memory
  }
  meta {
    author: "Dieter Best"
    email: "Dieter.Best@cdph.ca.gov"
    description: "CDC TB profiler"
  }
  parameter_meta {
    # Inputs:
    R1: "Forward reads in read pair"
    R2: "Reverse reads in read pair"
    SAMPL: "Sample name"
  }
}

workflow WfVarPipeline {
  input {
    File R1
    File R2
    File REF
    String SAMPL
    Int threads = 1
    String genome = "NC_000962"
  }
  call RunVarpipeline {
    input:
    R1 = R1,
    R2 = R2,
    REF = REF,
    SAMPL = SAMPL,
    threads = threads,
    genome = genome
  }
}
