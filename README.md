# ouroborUSCO
Iterative extraction of busco loci from short read datasets for phylo purposes


esentially just:

```
for iter in 1..N do
  align reads to $ref
  call variants w/ mpileup
  update $ref with $variants
  extract target regions from $ref, set to new $ref 
```

Hopefully after ~5-10 iterations this converges on each sample's true genotype
