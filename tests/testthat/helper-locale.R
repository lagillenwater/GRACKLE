## TRUE when this machine can collate in `locale`, so tests that compare collations can run.
collationAvailable <- function(locale = "en_US.UTF-8") {
  identical(suppressWarnings(withr::with_collate(locale, Sys.getlocale("LC_COLLATE"))), locale)
}
