# GitLab CI

This text documents the GitLab continuous integration (CI) YAML script (`.gitlab-ci.yml`) used in GVEC.

Currently, the script invokes a pool of runners that are shared among all MPCDF users. The novel module-enabled MPCDF Docker images infrastructure allowa to choose different Docker images.

In this document we analyse the source code of GVEC's CI script. The first subsection ([Organisation](#organisation)) provides an overview of the structure of the script, which will be useful to keep in mind throughout the document. The second subsection ([Jobs & templates per stage](#jobs-templates-per-stage)) provides more implementation details. We analyse the jobs and their templates together according to the stage they belong to, in order to make the description easier to follow. However, it should be noted that in the actual script all templates appear first, followed by the jobs, organised by stage. Finally, in the last subsection ([Conclusion](#conclusion)) a few basic concluding remarks are given. At the very end, a list of links to GitLab YAML reference pages is given, corresponding to concepts or keywords used throughout the text.

## Organisation

The CI script is organised into three main sections, as described here.

### Stage declaration & global variables

As mentioned in the introduction, the CI jobs are grouped into four different [`stages:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#stage). They are declared at the beginning of the CI script, and looks like this

```yaml
stages:
  - env
  - build
  - run
  - regression
  - post
```

There are also some global variables, assigned as follows

```yaml
# ========================================================================
# GLOBALS
# ========================================================================

variables:
  GIT_STRATEGY: none
  HASH_TAG: $CI_COMMIT_REF_NAME
  HASH_TAG_REFERENCE: develop
  PYTEST_EXEC_CMD: "python -m pytest -v -r A"
  PYTEST_MARKER_OPTS: "example"
  PYTEST_KEY_OPTS: ""
```

of which the last lines are related to the testing framework used in GVEC and the two variables before that allow to specify two different commits, the one triggering the CI pipeline and another one that in the example above corresponds to a tag commit. More about this topic later on, when we discuss the `regression` stage.

### Templates

An extensive use of CI job templates (named as `.tmpl_<name>` in the script) is made to reduce code repetition. These templates are configuration sections that can be reused by several CI jobs by invoking the [`extends:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#extends) keyword. Therefore they simplify the source code and improve maintainability, since future changes to the CI script are necessarily more localised. In terms of the organisation of the GVEC CI script, templates appear before any of the stages.

### Jobs

The next section of the CI script comprises the CI jobs specification within each CI stage. These jobs extend the templates declared before, and therefore inherit the corresponding configuration.

Below, we use (often reduced) CI script snippets to explain the different kinds of CI templates used by GVEC, and how they are extended and combined by the CI jobs within each CI stage.

## Jobs & templates per stage

The code snippets shown in the following paragraphs attempt to illustrate major aspects of GVEC's CI script. Some of the actual configuration details are skimmed down (wherever "(...)" appears) to maximise clarity. Additionally, all examples refer to the novel MPCDF Docker images, for simplicity. This is because the corresponding implementation for the remaining CI images/runners is quite similar. For the complete implementation details, the reader is advised to look into the CI script `.gitlab-ci.yml` directly.


### Templates independent of stage

The first two templates that appear in the CI script are independent of any of the script stages, in the sense that they are reused in several of them, as we will see later.

```yaml
# ========================================================================
# TEMPLATES INDEPENDENT OF STAGE
# ========================================================================


# choose the MPCDF Docker image intel_2023_1_0_x:latest on a shared runner
.tmpl_mpcdfci_intel2023:
  image: gitlab-registry.mpcdf.mpg.de/mpcdf/ci-module-image/intel_2023_1_0_x:latest
  tags:
    - shared
  variables:
    CMAKE_HOSTNAME: "mpcdfcirunner"
    CURR_CMP: "intel"


# setup software stack (load modules)
.tmpl_before_script_modules:
  before_script:
    - . ./CI_setup/${CMAKE_HOSTNAME}_setup_${CURR_CMP}
```

The first template (`.tmpl_mpcdfci_intel2023`) shows how to select a particular MPCDF Docker image, which an example provided corresponds to an MPCDF Docker image (`intel_2023_1_0_x:latest`).

The second template (`.tmpl_before_script_modules`) contains a  [`before_script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#before_script) section, which lists an overriding set of commands that are executed before a CI job. In this case, the only command listed sets up the software stack for a particular host, i.e. CI image/runner, specified by the variable `${CMAKE_HOSTNAME}_setup_${CURR_CMP}`. In particular, in this case it sources the following Shell script (`CI_setup/mpcdfcirunner_setup_intel`):

```bash
#!/bin/bash

module purge
module load git git-lfs cmake
module load intel mkl anaconda
module load hdf5-serial netcdf-serial
module list
export FC=`which mpiifort`
export CC=`which mpiicc`
export CXX=`which mpiicx`
```

### Templates for stage `env`

The next section of the CI script has a template for the `env` stage. It provides the Shell script to be executed by a runner, i.e. the [CI YAML `script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#script) keyword, which in this case comprises sending the values of system environment variable to stdout, for future reference.

```yaml
# ========================================================================
# TEMPLATES FOR STAGE "env" 
# ========================================================================

# env related shell commands to stdout
.tmpl_script_env:
  script:
    - echo "Pipeline environment for branch:" $CI_COMMIT_REF_NAME
    - printenv
    (...)
```

Jobs that extend any of the previous templates will inherit the specified configuration, as we will see next.

### Stage `env`

Here we skip the next set of templates for the subsequent stages and jump into directly the stage `env` to illustrate how the templates introduced so far are extended in practice. The corresponding job is

```yaml
# ========================================================================
# Stage "env"
# ========================================================================

# printout MPCDF Docker (Intel 2023) shared runner environment
mpcdfci_intel2023_env:
  stage: env
  extends:
    - .tmpl_mpcdfci_intel2023
    - .tmpl_script_env
```

which chooses the CI image/runner by inheriting from the template `.tmpl_mpcdfci_intel2023` and executes the commands specified in `.tmpl_script_env` that output environment variables values to stdout.

If, for instance, we would like output the environment variables values for a different CI image, we simply need to duplicate the job above, rename it and replace the template entry `.tmpl_mpcdfci_intel2023` under `extends:` with another one specifying the desired Docker image.


### Templates for stage `build`

Now we jump back to the template section of the CI script to describe the templates that are going to be extended in the the `build` stage.

```yaml
# ========================================================================
# TEMPLATES FOR STAGE "build" 
# ========================================================================

# GVEC CMake variables
.vars_cmake_def_opts:
  variables:
    COMPILE_GVEC: "ON"
    LINK_TO_NETCDF: "ON"
    COMPILE_CONVERTERS: "ON"
    (...)

# combination of build variable values
.vars_matrix_build:
  parallel:
    matrix:
      - CMP_MODE: ["Debug", "Release"]
        OMP_MODE: ["ompOFF", "ompON"]
  variables:
    MPI_MODE: "mpiOFF"

# target build directory defined by unique variable BUILDNAME
.tmpl_setup_build:
  variables:
    BUILDNAME: ${HASH_TAG}_${CURR_CMP}_${CMP_MODE}_${OMP_MODE}_${MPI_MODE}
    (...)

# preparation before build job
.tmpl_before_script_build:
  before_script:
  	(...)
    - if [ ${HASH_TAG} != ${CI_COMMIT_BRANCH} ]; then git checkout ${HASH_TAG}; fi
    - rm -rf build_${BUILDNAME}; mkdir -p build_${BUILDNAME}
    - cd build_${BUILDNAME}; pwd
	(...)

# build job script
.tmpl_script_build:
  script:
    - cmake ["CMAKE options"] ../.
	(...)
  artifacts:
    name: "${CI_PIPELINE_ID}_${BUILDNAME}"
    paths:
      - build_${BUILDNAME}
      (...)
```

The first template (`.vars_cmake_def_opts`) sets the default values for CMake variables needed to build GVEC.

The second template (`.vars_matrix_build`) sets up additional variables (`CMP_MODE` and `OMP_MODE`), but in a way to yield all possible combinations of their values using the [`parallel:matrix:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#parallelmatrix) keyword. This allows inheriting jobs to be executed by the CI runners concurrently, so it helps improving the CI pipeline performance.

The third template (`.tmpl_setup_build`) sets up variables that are relevant for this stage, including one that specifies a unique directory name for each CI job, where GVEC will be built.

The forth template (`.tmpl_before_script_build`) provides a [`before_script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#before_script) section that, as we have seen before, lists an overriding a set of commands that are executed before a CI job. These include in this case, creating the specified build directory, but also making an important distinction between using the GVEC version from the commit triggering the CI pipeline or another previous commit, corresponding to a tag or a release version of GVEC. This topic shall be discussed in more detail later on, when we describe the `regression` stage.

The last template (`.tmpl_script_build`) provides the [`script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#script) section for an inheriting CI job. It contains a list of Shell commands to be executed, which in this case build a configuration of GVEC specified by the variables defined in the first two templates. The latter specify things like choosing between building a debug or release version, with which compiler, etc. Finally, the last section in this template, [`artifacts:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#artifacts), specifies a list of files and directories that will be attached to the job upon success (by default). These, which here include the GVEC executable produced after a successful build job, can subsequently be used by dependent CI jobs, like the ones in the `run` stage.


### Stage `build`

Like we did before, here we jump forward to the corresponding `build` stage, in order to illustrate how a typical GVEC CI `build` job extends the previous templates. We use a CI job example that executes on the MPCDF shared runner using the same Docker image as before, namely:

```yaml
# ========================================================================
# Stage "build"
# ========================================================================

# build with MPCDF Docker (Intel 2023) shared runner environment
mpcdfci_intel2023_build:
  stage: build
  before_script:
    - !reference [".tmpl_before_script_modules", "before_script"]
    - !reference [".tmpl_before_script_build","before_script"]
  extends:
    - .tmpl_mpcdfci_intel2023
    - .tmpl_setup_build
    - .tmpl_script_build
    - .vars_cmake_def_opts
    - .vars_matrix_build
  needs: [mpcdfci_intel2023_env]
```

The first section is a [`before_script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#before_script). It allows the job to inherit from the very first template mentioned in this document, `.tmpl_before_script_modules`, which essentially loads the necessary modules. It also inherits from the template `.tmpl_before_script_build`, presented in the previous subsection, which executes commands in preparation of the build script, like creating the build target directory. Note that the keyword `!reference` is used instead of `extends:` due to the overriding property of the `before_script:` sections, otherwise it would not be possible to inherit from multiple of these sections without resetting some their variable values.

The next section is an [`extends:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#extends) allowing this CI job to inherit from five different templates, all covered before. They basically specify the desired Docker image (`.tmpl_mpcdfci_intel2023`), define the different target build directories (`.tmpl_setup_build`) and execute build commands (`.tmpl_script_build`) for all versions of GVEC specified by the variables inherited from the templates `.vars_cmake_def_opts` and `.vars_matrix_build` together with the ones specified in the job section [`variables:`]{https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#variables}. Recall that, because the template `.vars_matrix_build` uses the keyword [`parallel:matrix:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#parallelmatrix) together with two variables each having two possible values, the build job `mpcdfci_intel2023_build` generates in practice four different CI jobs that can be executed concurrently in the pipeline, which improves performance.

Finally, the last section specifies the dependency of these jobs on the previous job `mpcdfci_intel2023_env `. By using here the keyword [`needs:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#needs) creates a [direct acyclic graph](https://docs.gitlab.com/ee/ci/directed_acyclic_graph/) between the different jobs, meaning that the different stages (`env`, `build`, `run` and `regression`) do not have to be executed sequentially. This naturally improves also CI pipeline performance.


### Templates for stage `run`

After each job from the previous stage (`build`) is finished, the corresponding job(s) in the stage `run` can start. As always, at this point we step back a little and analyse first the templates that are specific to the stage under consideration.

```yaml
# ========================================================================
# TEMPLATES FOR STAGE "run"
# ========================================================================

# combination of runtime variable values
.vars_matrix_run:
  parallel:
    matrix:
      - CMP_MODE: ["Debug"]
        OMP_MODE: ["ompOFF","ompON"]
      - CMP_MODE: ["Release"]
        OMP_MODE: ["ompON"]
    (...)

# target run directory defined by unique variable CASENAME
.tmpl_setup_run:
  variables:
    BUILDNAME: ${HASH_TAG}_${CURR_CMP}_${CMP_MODE}_${OMP_MODE}_${MPI_MODE}
    CASENAME: ${HASH_TAG}_${CURR_CMP}_${CMP_MODE}_${OMP_MODE}_${MPI_MODE}${MPI_RNKS_MODE}
    (...)

# preparation before run job
.tmpl_before_script_exportvars:
  before_script:
    - echo "BUILDNAME is ${BUILDNAME} >> build_${BUILDNAME}/"
    - echo "CASENAME is ${CASENAME} >> CIrun_${CASENAME}/"
    (...)

# more preparation before run job
.tmpl_before_script_run:
  before_script:
    - export PYTEST_DIR_OPTS="--builddir=build_${BUILDNAME} --rundir=CIrun_${CASENAME}"
    - rm -rf CIrun_${CASENAME}
	(...)

# run job script
.tmpl_script_run:
  script:
    - ${PYTEST_EXEC_CMD} --log-file=log_pytest_run_norestart.txt -m "${PYTEST_MARKER_OPTS} and run_stage and (not restart)" $PYTEST_DIR_OPTS -k "${PYTEST_KEY_OPTS}"
  artifacts:
    name: "${CASENAME}"
    paths:
      - CIrun_${CASENAME}
	(...)
```

Since these templates are very similar to the ones from the stage `build`, we are not going to repeat the description made therein, but rather highlight their differences. In particular, we note that there is an extra variable in `.tmpl_setup_run`, namely, `CASENAME`. This specifies the target directory where specified version of GVEC is going to run, which is different from the target build directory `BUILDNAME`. Naturally, the latter is also needed, because that is were the GVEC executable compiled during the `build` stage is stored. Another clear difference is, of course, the commands executed in the `script:` section, since now job runs GVEC, instead of building it. Moreover, looking more carefully at the last Shell command of the `.tmpl_script_run` template, GVEC is not executed directly. Rather, a Python framework based on Pytest is invoked, which then executes several GVEC cases by providing the corresponding input parameters, and further compares output results to stored reference values. In other words, it performs already a set of end-to-end tests. CI Artifacts are also used use to save the results for the subsequent stages (`regression` and `postprocessing`).


### Stage `run`

The corresponding stage `run` extends these templates, pretty much the same way the stage `build` extended its counterpart templates. Furthermore, the `run` stage CI job specified below is very similar to the one we discussed in the `build` stage previously. So, once more, we'll focus only on their differences, which in this case, is essentially just one, namely, the [needs:](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#needs) keyword. It additionally specifies which [`artifacts:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#artifacts) are required. Notably in this case, the whole target build directory produced during job `mpcdfci_intel2023_build ` of the previous `build` stage, which includes the corresponding GVEC executable.

```yaml
# ========================================================================
# Stage "run"
# ========================================================================

# run with MPCDF Docker (intel2023) shared runner environment
mpcdfci_intel2023_run:
  stage: run
  before_script:
    - !reference [".tmpl_before_script_modules", "before_script"]
    - !reference [".tmpl_before_script_exportvars", "before_script"]
    - !reference [".tmpl_before_script_run", "before_script"]
  extends:
    - .tmpl_mpcdfci_intel2023
    - .tmpl_setup_run
    - .tmpl_script_run
    - .vars_matrix_run
  needs:
    - job: mpcdfci_intel2023_build
      artifacts: true
```


### Templates for stage `regression`

We have reached the last stage in the GVEC CI script, namely, `regression`. Why do we have this additional stage, when regression tests were already performed in the `run` stage? The answer is that we want to be able to compare output of different GVEC builds, either from the same commit (e.g. MPI vs. no-MPI, compilation in debug vs. release mode, or compilation with Intel vs. GNU compilers), or between different commits (e.g. current branch vs. tag commit). In the `run` stage we could only compare current results with reference stored values, which implies certain limitations. Most notably, there is no control over the conditions under which the runs producing the reference data were made. Conversely, when comparing two different commits of GVEC, compiled on-the-fly with the same software stack and executed on the same hardware, one is sure that any discrepancies that arise in the test results are due to source code differences between them, and not e.g., because of different compiler vendors or optimisation flags were used. Alternatively, being able to compare the same GVEC commit built with different compiler vendors or with built the same compiler but using different options (e.g. OpenMP vs hybrid OpenMP/MPI), gives information about the impact of those choices. There are of course many other possible comparisons that can be made, and overall, the idea is that, being able to make all these on-the-fly makes the regression tests much more powerful.

Below we are going to analyse how the corresponding CI stage is implemented in practice, starting first, as we have done for the previous stages, to inspect the corresponding templates. Incidentely, there is only one template that is specific to the `regression` stage, namely

```yaml
# ========================================================================
# TEMPLATES FOR STAGE "regression"
# ========================================================================

# regression testing job script
.tmpl_script_regression:
  script:
    - ${PYTEST_EXEC_CMD} --log-file=log_pytest_regression.txt -m "${PYTEST_MARKER_OPTS} and regression_stage" $PYTEST_DIR_OPTS
  artifacts:
    name: "${CASENAME_1}_vs_${CASENAME_2}"
    paths:
      - CIrun_${CASENAME_1}
      - CIrun_${CASENAME_2}
      (...)
```

This contains mostly the command to execute a Pyhton wrapper in the [`script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#script) section, and declares which files are to be stored as artifacts, for future reference.


### Stage `regression`

We now show an example of couple of `regression` jobs extending the template above, together with other templates that were presented previously. At first glance, ones sees that there are two sets of variables with suffixes `_1` and `_2`. These correspond to each of the two `run` stage jobs that are compared.

```yaml
# ========================================================================
# Stage "regression"
# ========================================================================

# compare results between intel and gnu on current branches (no MPI)
mpcdfci_intel-vs-gnu_reg:
  stage: regression
  extends:
    - .tmpl_mpcdfci_intel2023     # need to choose one (does not matter which) Docker image
    - .tmpl_before_script_modules # to be able to load modules
    - .tmpl_script_regression
    - .vars_matrix_run
  variables:
    HASH_TAG_1: ${HASH_TAG}
    HASH_TAG_2: ${HASH_TAG}
    CURR_CMP_1: intel
    CURR_CMP_2: gnu
    CASENAME_1: ${HASH_TAG_1}_${CURR_CMP_1}_${CMP_MODE}_${OMP_MODE}_mpiOFF
    CASENAME_2: ${HASH_TAG_2}_${CURR_CMP_2}_${CMP_MODE}_${OMP_MODE}_mpiOFF
  needs:
    - job: mpcdfci_intel2023_run
      artifacts: true
    - job: mpcdfci_gcc13_run
      artifacts: true

# compare results between current and tag branches (no MPI)
mpcdfci_intel-vs-tag_reg:
  stage: regression
  extends:
    - .tmpl_mpcdfci_intel2023     # need to choose one (does not matter which) Docker image
    - .tmpl_before_script_modules # to be able to load modules
    - .tmpl_script_regression
    - .vars_matrix_run
  variables:
    HASH_TAG_1: ${HASH_TAG}
    HASH_TAG_2: ${HASH_TAG_REFERENCE}
    CURR_CMP_1: ${CURR_CMP}
    CURR_CMP_2: ${CURR_CMP}
    CASENAME_1: ${HASH_TAG_1}_${CURR_CMP_1}_${CMP_MODE}_${OMP_MODE}_mpiOFF
    CASENAME_2: ${HASH_TAG_2}_${CURR_CMP_2}_${CMP_MODE}_${OMP_MODE}_mpiOFF
  needs:
    - job: mpcdfci_intel2023_run
      artifacts: true
    - job: mpcdfci_intel2023_run_tag
      artifacts: true
```

The first job (`mpcdfci_intel-vs-gnu_reg`) compares the same commit of GVEC when built using Intel and GNU compilers (no MPI). The two first templates extended (`.tmpl_mpcdfci_intel2023 ` and `.tmpl_before_script_modules `) are from previous stages. As we know already, they choose a Docker image and load the respective modules. Note that at this level, which image is chosen is not so important. What matters is that the module for `anaconda` is loaded, since Pyhton >3.10 is needed for the Pytest framework, which is called from the template `tmpl_script_regression`. This wrapper compares output files from the two different directories specified by the variables `CASENAME_1` and `CASENAME_2`, which are in turn made by two sets of values of build variables: `CURR_CMP`, `HASH_TAG`, `CMP_MODE` and `OMP_MODE`. The latter two are declared via the [`parallel:matrix:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#parallelmatrix) keyword via the template `.vars_matrix_run` which triggers three different jobs in parallel. Naturally, this must be consistent with the set of jobs executed during the `run` stage, since now we depend on the existence of results produced there. In other words, the combination of variables specified here (`regression` stage) must be contained in the respective combination of variables specified at the `run` stage. This is guaranteed here by using the same parallel matrix template in both stages. Finally, we need to explicitly prescribe a dependency on `run` stage jobs, which must be chosen consistently with the values picked for those same variables. In this case, those jobs are `mpcdfci_intel2023_run` and `mpcdfci_gcc13_run`, together with their artifacts, which include the output files to be compared.

The second set of `regression` jobs shown (`mpcdfci_intel-vs-tag_reg`) is very similar, but with some notable differences. Namely, the values for all `_1` and `_2` variables are the same, except for `HASH_TAG_1` and `HASH_TAG_2`, which therefore specify different commits. This means that, instead of comparing the results of the same commit built with two compilers, like before, we now compare two different commits using the same build setup. In particular, `HASH_TAG_1` refers to the commit triggering the pipeline (`HASH_TAG`), and the `HASH_TAG_2` to the specific commit tag `HASH_TAG_REFERENCE=v0.2.0`. The latter was specified as a global variable, as mentioned at the beginning of this document. In order to be consistent with this specification of a different commit, the `run` job dependencies and corresponding artifacts now include the job `intel2023_run_tag` from the stage `run`, which itself depends on a `build` stage job that checks out a commit specified by the variable `HASH_TAG_REFERENCE`. Therefore, it is now time to revisit the corresponding Shell command, which is part of the `.tmpl_before_script_build` template, presented when we discussed the stage `build` and its templates. Namely,

```yaml
.tmpl_before_script_build:
  before_script:
    (...)
    - if [ ${HASH_TAG} != ${CI_COMMIT_REF_NAME} ]; then git fetch --tags; git lfs fetch; git checkout ${HASH_TAG}; fi
    (...)
```

Here we see that, if the variable `HASH_TAG`, which is globally set to be `CI_COMMIT_REF_NAME` by default (as we saw before), is subsequently set to something else, then an explicit `git checkout` of the corresponding commit is issued. This commit is then used for the remaining CI jobs, which is what enables regression tests between different commits in GVEC's CI pipeline.

Finally, it is noteworthy to mention that the global variable `HASH_TAG_REFERENCE` can be specified externally to override the value hardcoded in the CI script (in this case the branch `develop`). E.g., this can be achieved by creating [scheduled pipelines](https://gitlab.mpcdf.mpg.de/help/ci/pipelines/schedules) that are triggered automatically, e.g. at regular intervals, and specifying this variable and its desired value in its configuration.


### Templates for stage `postprocessing`


### Stage `postprocessing`


## Conclusion

As we have seen, all the jobs within a given CI stage are very similar to each other, in part thanks to the usage of common templates. In this document, whose purpose is to describe the machinery behind GVEC's CI script, only a subset of the CI jobs implemented are covered, and moreover, in a simplified fashion. For the complete implementation details, including all additionally available cases that were not shown here, the reader is advised to inspect directly the script file `.gitlab-ci.yml`.


## References

##### List of references to keywords made throughout the document:

* [`stage:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#stage)
* [`extends:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#extends)
* [`before_script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#before_script)
* [`script:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#script)
* [`parallel:matrix:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#parallelmatrix)
* [`artifacts:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#artifacts)
* [`needs:`](https://gitlab.mpcdf.mpg.de/help/ci/yaml/index#needs)

##### List of reference to Git related concepts or services mention in this document:

* [direct acyclic graph](https://docs.gitlab.com/ee/ci/directed_acyclic_graph/)
* [scheduled pipelines](https://gitlab.mpcdf.mpg.de/help/ci/pipelines/schedules)
