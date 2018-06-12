# Using the Mammouth HPC server

This document provides information on how to get started with using the Mammouth HPC.

## Connecting to the Mammouth HPC scheduler:

Host: workshop2018a.vhost37.genap.ca

Login / Password: Provided to you on a separate piece of paper at the beginning of the workshop.

### Linux/Mac

Open a terminal and use the following command from the command line: ssh studXX@workshop2018a.vhost37.genap.ca The username and password are the ones that we provided on a separate piece of paper.

### Windows

Connect to the studXX@workshop2018a.vhost37.genap.ca server using the Putty tool that you installed for this workshop.

## Once you are logged in

When you are logged on a Compute  Canada HPC server you are on a logging node which is only dedicated to data transfer, to run very small jobs and to launch job on the scheduler.  

In this workshop environment you are directly working in a working node. So you will not have to use the scheduler system. 

### Loading a module

All available software are organized as modules that you can load. For example, to start using Java 8, you need to load the Java 8 module.

1- Get the list of available modules  

`[lect01@workshop103 ~]$ module avail`

2- Load your module of choice  

`[lect01@workshop103 ~]$ module load mugqic/java/openjdk-jdk1.8.0_72`

3- You can now run the module you loaded directly  

`[lect01@workshop103 ~]$ javac -version `

> javac 1.8.0_72-ea


## Downloading a file from Mamouth server to your computer

### Linux/Mac

Use the following from the command line on your local computer:  

```
cd /path/to/destination/folder 
scp studXX@workshop2018a.vhost37.genap.ca:path/to/destination/myfile.txt .
```


To download a full directory, add the -r option 

```
cd /path/to/destination/folder 
scp -r studXX@workshop2018a.vhost37.genap.ca:path/to/destination/myfolder .
```


### Windows

Use the WinSCP software that you installed for this workshop, connecting to`workshop2018a.vhost37.genap.ca` .




## (Optional) Launching a job on the scheduler

Note: Complete documentation on the scheduler is available [here](https://wiki.calculquebec.ca/w/Ex%C3%A9cuter_une_t%C3%A2che/en#tab=tab7)

Jobs are submitted to the scheduler using the qsub command. The important qsub parameters are the following: 

 * nodes : total number of machines to use for this job 
 * ppn : total number of cores to use for this job 
 * -d workdir : The working directory, in other words the path in which this command will be executed

Please keep in mind that you always need to load the modules you want to use in your scripts.

### Using a bash script

You can specify a bash script to launch as a compute job directly from the command line, like the `testjob.sh` script in this example: `qsub -l nodes=1:ppn=1 -A bem-651-ae testjob.sh`

where the content of testjob.sh could be: 

``` 
#!/bin/bash 
# PBS -l nodes=1:ppn=2

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 

java -Xmx8G -jar ${BVATOOLS_JAR} readsqc --read1 raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz --read2 raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz --threads 2 --regionName ACTL8 --output originalQC/

```

### Specifying the whole command from the prompt

If you do not want to put your command in a script, you can write it directly in the prompt using the echo command like this: 

```
echo "touch test.txt" | qsub -l nodes=1:ppn=1 -d .
```

### Viewing the current state of your jobs

Use the following command to see the current state of the jobs you launched: showq -u <LOGIN_NAME> -n -v

 * Idle - means the job is waiting for its turn in the queue
 * Running - means the job is currently running on one of the execution nodes

It can sometimes take up to a minute for your jobs to appear in the showq output.
