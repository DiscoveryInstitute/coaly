# Coaly
Coaly is a simple simulation tool for doing quick-and-dirty calculations of expected Allele Frequency Spectra to test hypotheses quickly under different population genetics scenarios. There are three methods included, which are described in the publication: [insert doi reference here].

### COMPILE, SET-UP, AND RUN
---
**Compile** using `./build.sh` to create an executable `Coaly`.
(Make sure you have the latest version of g++).
Copy the executable to wherever you will run it. 

**Edit** input files. The first file, which you could name `demography.txt` contains the haploid population size over time.
The second file, suggested name `mutationrate.txt` contains the mutation rate over time. 
The third file, suggested name `primordial.txt`, is optional. It describeds primordial variation, if present.
The format of these three files is described below. 

**Run** with `./Coaly method demography.txt mutation.txt [primordial.txt] > output.txt`.
The square brackets `[]` around the fourth argument to the program indicate that it is optional.
The parameter `method` can be one of
* `stochastic` - the stochastic simulation method described in the reference.
* `forward` - the accurate forward matrix calculation method described in the reference.
* `backward` - the fast approximate backward coalescent method "Coaly" described in the reference.


The file names ending `.txt` are only suggested names. 

### DEMOGRAPHIC-HISTORY FILE
---
This file contains a series of pairs: a generation number, and then a haploid population size. 
The size of the population is interpolated between each point. 
The generation numbers are counted back in time, but they are listed forwards in time from the founding generation to the present. Therefore the generation numbers must be decreasing. For example, for a simulation of 5000 generations with a constant haploid population of 1000 :
```
5000 1000 
0    1000 
```
For a population that slowly grows from 1000 to 2000 and then plateaus at 2000:
``` 
5000 1000
4000 2000
0    2000
```

### MUTATION-RATE-HISTORY FILE
---
This file is very similar to the demographic-history file except that it has a series of pairs: a generation number and the mutation rate at that generation. Mutation rates are again linearly interpolated between points. 
For example, for mutation rate that falls linearly from 2 x 10^-8 to 1 x 10^-8 over 2000 generations:
```
2000 2e-8
0 1e-8
```

### PRIMORDIAL-VARIATON FILE
---
This file has a different format from the previous two. It also contains pairs: a frequency, and then the number of primordial variants having that frequency. The frequency is an integer: the total number of members of the haploid population that have the variant. It must not exceed the total founding population. 
For example, for a diploid founding couple, the haploid population is size 4, and the primordial variation may be composed of variants of frequency 1 or 2:
```
1 1e6
2 5e5
```

---
---
