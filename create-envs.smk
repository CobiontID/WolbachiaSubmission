my_envs = ['busco.yaml',
           'cdhit.yaml',
           'hifiasm.yaml',
           'hmmer.yaml',
           'prokka.yaml',
           'ete3.yaml',
           'treebuild.yaml']

rule make_all_envs:
    input:
        expand("created-{name}", name=my_envs)

for env_file in my_envs:
    rule:
        output:
            temp("created-%s" % env_file)
        conda:
            "envs/%s" % env_file
        shell:
            "touch {output}"