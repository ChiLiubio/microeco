
.onAttach <- function(...) {
	withr::with_preserve_seed({
		# Show messages only in interactive sessions to prevent interference during non-interactive checks (e.g., R CMD check)
		if (!interactive() || stats::runif(1) > 0.1) {
			return()
		}
		
		tips <- c(
			"For microeco version 2, please cite: https://doi.org/10.1002/imt2.70132", 
			"Run ?microeco to find links to GitHub Issues and Tutorials"
		)
		
		tip <- sample(tips, 1)
		packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
	})
}
