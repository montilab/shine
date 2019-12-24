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
  path.eset = "{path_eset}"
  path.genes = "{path_genes}"
  path.blanket = "{path_blanket}"
  outdir = "."
  iters = "{iters}"
  cores = "{cores}"
  condition = "{condition}"
  email = ""
}}
'''.format(**kwrgs))

def nf_head():
    return('''\
#!/usr/bin/env nextflow

log.info """\
-
P I P E L I N E ~ Configuration
===============================
eset      : ${params.path.eset}
genes     : ${params.path.genes}
blanket   : ${params.path.blanket}
outdir    : ${params.outdir}
iters     : ${params.iters}
cores     : ${params.cores}
condition : ${params.condition}
-
"""

eset = file(params.path.eset)
genes = file(params.path.genes)
blanket = file(params.path.blanket)
''')

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

def nf_process(label, include, blanket=False, prior=False):
    process = """
process {0} {{
    cache "deep"
    publishDir "$params.outdir/networks", pattern: "*.rds", mode: "copy"
    publishDir "$params.outdir/logs", pattern: "*.log", mode: "copy"

    input:
    file eset
    file genes
    val include from "{1}"
    val label from "{0}"
    {2}
    {3}

    output:
    file '*.rds' into {0}_rds
    file '*.log'

    script:
    template 'run.sh'
}}
"""
    if blanket: blanket = 'file blanket'
    else:       blanket = 'val blanket from "/"'

    if prior: prior = 'file prior from {0}_rds'.format(prior)
    else:     prior = 'val prior from "{0}"'.format("/")

    return( process.format(label, include, blanket, prior) )

def nf_workflow(hierarchy, blanket):

    workflow = [nf_head()]

    leaves = {}
    for i in hierarchy.split('\n'):
        if i == "":
            continue 

        # Networks outside of a hierarchy
        if '->' not in i:
            label = i.strip(' ')
            include = split_label(label)
            workflow.append( nf_process(label, include, blanket) )
            continue

        # Networks in a hierarchy
        s = i.split('->')
        leaves[s[1].strip(' ')] = s[0].strip(' ')   

    # Roots first
    tree = set([i for i in leaves.values() if i not in leaves.keys()])
    for label in tree:
        include = split_label(label)
        workflow.append( nf_process(label, include, blanket) )

    while len(leaves) > 0:
        for child, parent in leaves.items():
            if parent in tree:
                include = split_label(child)
                workflow.append( nf_process(child, include, blanket, parent) )
                tree.add(child)
        
        for i in tree: 
            try: del leaves[i]
            except: pass

    workflow.append(nf_tail())
    
    return(workflow)

def nf_build(outdir,
             hierarchy,
             condition,
             path_eset,
             path_genes,
             path_blanket,
             iters,
             cores):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if path_blanket == "/":
        blanket = False
    else:
        blanket = True

    with open(os.path.join(outdir, "hierarchy.nf"), "w") as outfile:
        for i in nf_workflow(hierarchy, blanket):
            outfile.write(i)

    with open(os.path.join(outdir, "hierarchy.config"), "w") as outfile:
        outfile.write( nf_config(condition=condition,
                                 path_eset=path_eset,
                                 path_genes=path_genes,
                                 path_blanket=path_blanket,
                                 iters=int(iters),
                                 cores=int(cores)) )
