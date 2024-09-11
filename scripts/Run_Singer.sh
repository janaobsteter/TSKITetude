vcfDir=
vcf=
ratio=1
mutation_rate=
output=
num_samples=10
thinning=10
num_cores=
path_to_singer=
#Ne=



${path_to_singer}/parallel_singer  -m ${mutation_rate} \ #-Ne ${Ne}
-vcf ${vcfDir}/${vcf} -output ${output} \
-n ${num_samples} -thin ${thinning} -ratio ${ratio} -num_cores ${num_cores} 