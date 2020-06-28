# This function is used to translate AM probe IDs to OGS 3.2 genes using the AM_to_OGS3_2.txt table

translate <- 
    function(to_translate, translation_table)
    {
        to_translate %>%
        match(translation_table[,1],nomatch=NA) %>%
        sapply(
            function(x)
            {
                if (!is.na(x))
                    translation_table[x,2]
                else 
                    NA
            }
        )
    }