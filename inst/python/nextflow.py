import os

def split_label(x):
    return(" ".join([i for i in x.split('_')]))

def nf_config(**kwrgs):
    return('''\
profiles {{
  local {{includeConfig 'configs/local.config'}}
  docker {{includeConfig 'configs/docker.config'}}
  sge {{includeConfig 'configs/sge.config'}}
  aws {{includeConfig 'configs/aws.config'}}
}}

params {{
  data = "$baseDir/data/data.rds"
  output = "$baseDir"
  email = ""
}}
'''.format(**kwrgs))

def nf_head(**kwrgs):
    return('''\
#!/usr/bin/env nextflow

log.info """\
-
W O R K F L O W ~ Configuration
===============================
data      : ${{params.data}}
output    : ${{params.output}}
-------------------------------
Hierarchy
{hierarchy}
-
"""

data = file(params.data)
run = file("$baseDir/scripts/run.R")

'''.format(**kwrgs))

def nf_tail():
    return('''
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
  def subject = 'Estimation of Network Hierarchy Complete'
  def recipient = (params.email)
    ['mail', '-s', subject, recipient].execute() << """
    Pipeline Summary
    ---------------------------
    Timestamp    : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit Status  : ${workflow.exitStatus}
    Error Report : ${workflow.errorReport ?: '-'}
    """
}
''')

def nf_process(label, include, prior=False):
    process = """
process {0} {{
    cache "deep"
    publishDir "$params.output/networks", pattern: "*.rds", mode: "copy"
    publishDir "$params.output/logs", pattern: "*.log", mode: "copy"

    input:
    file data
    file run
    val include from "{1}"
    val label from "{0}"
    {2}

    output:
    file '*.rds' into {0}_rds
    file '*.log'

    script:
    template 'run.sh'
}}
"""

    if prior: prior = 'file prior from {0}_rds'.format(prior)
    else:     prior = 'val prior from "{0}"'.format("/")

    return( process.format(label, include, prior) )

def nf_workflow(hierarchy):

    workflow = [nf_head(hierarchy=hierarchy)]

    leaves = {}
    for i in hierarchy.split('\n'):
        if i == "":
            continue 

        # Networks outside of a hierarchy
        if '->' not in i:
            label = i.strip(' ')
            include = split_label(label)
            workflow.append( nf_process(label, include) )
            continue

        # Networks in a hierarchy
        s = i.split('->')
        leaves[s[1].strip(' ')] = s[0].strip(' ')   

    # Roots first
    tree = set([i for i in leaves.values() if i not in leaves.keys()])
    for label in tree:
        include = split_label(label)
        workflow.append( nf_process(label, include) )

    while len(leaves) > 0:
        for child, parent in leaves.items():
            if parent in tree:
                include = split_label(child)
                workflow.append( nf_process(child, include, parent) )
                tree.add(child)
        
        for i in tree: 
            try: del leaves[i]
            except: pass

    workflow.append(nf_tail())
    
    return(workflow)

def nf_build(hierarchy, data, outdir):

    with open(os.path.join(outdir, "workflow.nf"), "w") as outfile:
        for i in nf_workflow(hierarchy):
            outfile.write(i)

    with open(os.path.join(outdir, "workflow.config"), "w") as outfile:
        outfile.write(nf_config())
