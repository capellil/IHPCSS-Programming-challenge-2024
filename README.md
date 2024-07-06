# IHPCSS Programming challenge 2024

## Presentation slides
Tips and guidelines for the programming challenge are given in the presentation slides, available on the Moodle page at https://www.hpc-training.org/moodle/course/view.php?id=93#section-6.

## Run on compute nodes
You do not have GPUs on the login node. You **must** submit your job to compute nodes.

### Submit a job
If you do not use GPUs yet, you can submit your job to CPU nodes using the submission script: `submission_cpu_node.slurm`. However, if you use GPU programming, you must use the `submission_gpu_node.slurm` submission script to submit to a GPU node.

To submit a script, you need to type `sbatch your_submission_script.slurm`.

### Check for your job completion
To check the status of your current job(s), type `squeue -u $USER`. As a convenience, you can type `watch squeue -u $USER`: it will show the current list of jobs in queue and/or being executed, and refresh this list every 2 seconds by default.

## Graph generators
The graph processed is generated from inside the application. You will see that the code has two graph generators:
- `generate_nice_graph`: generates a graph that will easily offer parallelisation gains.
- `generate_sneaky_graph`: generates a graph that will require some analysis to get parallelisation gains. **This is the graph that will be used in submissions.**

## Submission
You are absolutely welcome to participate to the programming challenge for fun, and at your own pace. If you want to submit and your code to be taken into account for the competition, there are a few things to check to ensure fairness.

### Deadline
All submissions must be received at **Thursday 23:59 at the latest**. In other words, email received at Thursday 23:59 is fine, but Friday 00:00 is not.

### Teams
Teams must be made of between 1 and 3 individuals.

### Process
To submit, just send me an email at l.capelli@epcc.ed.ac.uk containing:
- your entire folder zipped (containing the submission script, makefile, everything)
- the list of team members.

## Questions?

- Is optimisation X allowed?
- Can I use more nodes?
- Can I change the compiler, the MPI process count, the OpenMP thread count?
- ...

Every time you wonder if something is allowed, feel free to ask. You are welcome to ask questions anytime. You can contact me by slack, by email at l.capelli@epcc.ed.ac.uk, or in person :)

**Remember, if you are doing the programming challenge just for fun and do not intend to submit at the end, you have no restriction and are free to tweak it as you please! :)**

## Stay tuned!

From experience, the repository regularly gets updates to clarify further which optimisations are and are not allowed, and having patches to fix any bugs that had done been spotted yet.

Therefore, you are **strongly encouraged** to regularly check this repository for updates, issuing `git pull` for instance.
Updates will also be posted on the slack channel https://ihpcss.slack.com/archives/C070MQHLPBQ.
