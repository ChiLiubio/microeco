
.onAttach <- function(...) {
	withr::with_preserve_seed({
		# Show messages only in interactive sessions to prevent interference during non-interactive checks (e.g., R CMD check)
		if (!interactive()) {
			return()
		}
		
		tips <- c(
			"MecoToolkit app for interactive microbiome analysis https://mecotoolkit.com/", 
			"For microeco v2, please cite: https://doi.org/10.1002/imt2.70132"
		)
		
		tip <- tips[1]
		packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
	})
}
