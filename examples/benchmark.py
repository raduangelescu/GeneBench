from genebench.benchmark import BenchmarkDiffMethods

benchmark = BenchmarkDiffMethods("config.json")
# run all methods and store results
benchmark.generate_method_results()
# generate the actual metrics and plots for comparison
benchmark.generate_comparisons()