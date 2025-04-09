
to_r_vec <- function(str) {
    str <- str_replace(str, "\\[", "\\(")
    str <- str_replace(str, "\\]", "\\)")
    str <- str_replace_all(str, "([^,\\(\\)]+)", "'\\1'")
    str <- paste0("c", str)
    out <- eval(parse(text = str))
    out <- str_trim(out)
    out
}
