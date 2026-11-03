import os, sys

template = open('ragtag.sh','r').read()
template_bigmem = open('ragtag_bigmem.sh','r').read()
count = 0
todo = open('../cut.txt','r').read().splitlines()

for fname in os.listdir('/data/blangme2/fasta/hprc-v2/final_update/fastas'):
    if fname.startswith('2'): # only unscaffolded assemblies
        f = fname.split('_')[1]
        idx, hap, _ = f.split('.')
        hap = '1' if hap == 'pat' or hap == 'hap1' else '2'
        name = idx + '.' + hap
        if name in todo:
        #if count <= 300:
        #    sh = template.format(name, fname)
        #else:
            sh = template_bigmem.format(name, fname)
            os.makedirs('./' + name, exist_ok=True)
            with open(name+'/'+name+'.sh','w') as out:
                out.write(sh)
        
            os.system(f'cd ./{name}; sbatch {name}.sh')
            count += 1



