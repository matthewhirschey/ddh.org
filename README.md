# ddh.com
data-driven hypothesis + drugs

## Data Generation
Generating data requires a Slurm cluster with singularity support.
To generate the data run:
```
./submit-slurm-jobs.sh
```

## Testing

### Running Tests
The tests can be run from the command line like so:
```
Rscript code/run_tests.R
```

### Single Module Test App
When working with a module running it in isolation may be preferable.
To help with this there is shiny app at [tests/TestShinyApp/app.R](tests/TestShinyApp/app.R).
To setup a page detail module for testing update:
- the data function
- the module render function
- the module server function

### Profiling Shiny Modules
To profile modules run the following from the command line:
```
Rscript tests/profile/profile.R
```
At the end the script will print out slow modules, files, and memory used.
