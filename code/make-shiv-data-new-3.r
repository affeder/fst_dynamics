#make shiv data


nucdat <- tbl_df(read.table("../dat/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE))
monkinf <- tbl_df(read.table("../dat/seqinfo.txt", header = TRUE, stringsAsFactors = FALSE))
aadat <- tbl_df(read.table("../dat/aminoacids.txt", header = TRUE, stringsAsFactors = FALSE))


dat <- bind_cols(monkinf, nucdat, aadat)
dat <- dat %>% filter(monk.id == "T98133" & (loc == "PLASMA" | loc == "LNRNA" | loc == "GUTRNA")) %>% 
               group_by(loc, sampweek) %>% mutate(index = 1:n())


#A site will be included in the variation if
# at least 90% of its reads are not "-"
# it has a minor allele of >1% of all data
# it appears at or before the first sampling timepoint with drug resistance
polymorphic <- function(x){

    if((x %>% summarize(tooManyDashes = mean(ident == "-") > 0.9))$tooManyDashes){ 
        minor = NA 
        highestfreq = NA
    } else{
        idents <- table((x %>%  filter(ident != "-"))$ident)
        identtab <- sort(idents, decreasing = TRUE)
        highestfreq <- names(identtab[1])
        minor <- sum(identtab[-1])
    }
    return(x %>% select(monk.id, pos) %>% distinct() %>% mutate(minor =  minor, hf = highestfreq))
    
}


x <- dat %>% filter(sampweek < 20) %>% filter(grepl("I|V", AA184)) %>%
    select(-starts_with('aa')) %>% 
    gather(pos, ident, -loc, -monk.id, -sampweek, -p.id, -f.id, -index) %>%
    group_by(pos) %>% do(polymorphic(.))


foreach(maf = c(1, 2, 5))%do%{

nucpos <- as.numeric(gsub("nuc", "", (x %>% filter(!is.na(minor), minor >= maf) )$pos))
aapos <- floor((nucpos - 1)/3)

tmpdat <- dat %>%  gather(pos, ident, -loc, -monk.id, -sampweek, -p.id, -f.id, -index)

#I'm going to just run this loopwise
dat.with.muts <- foreach(i = 1:length(nucpos), .combine = "bind_rows")%do%{

    nuccode <- paste0("nuc", nucpos[i])
    aacode <- paste0("AA", aapos[i])

    long <- dat %>% group_by(loc, sampweek, index) %>% 
        select(nuccode, aacode, index) %>% 
        gather(pos, ident, -loc, -sampweek, -index)

    refnuc <- names(sort(table((long %>% 
                   filter(sampweek <15, grepl("nuc", pos)))$ident), 
                   decreasing = TRUE)[1])
    refaa <- names(sort(table((long %>% 
                   filter(sampweek <15,grepl("AA", pos)))$ident), 
                   decreasing = TRUE)[1])
    
    long %>% mutate(ref = ifelse(grepl("nuc", pos), refnuc, refaa)) %>% 
        mutate(pos = gsub('[0-9]+', "", pos)) %>%
        mutate(nucmuts = ifelse((ident != ref & pos == "nuc") | (pos == "AA"), 
                                paste0(ref, aapos[i], ident), NA)) %>%
        select(-ident, -ref) %>% spread(pos, nucmuts) %>% 
        mutate(muts = ifelse(!is.na(nuc), AA, NA))

}


dattoprint <- dat.with.muts %>% 
    summarize(mutlist = paste0(muts, collapse = "-")) %>%
    mutate(mutlist = gsub("(NA-)+(NA)?", "", mutlist)) %>%
    mutate(mutlist = gsub("-$|(-NA)", "", mutlist)) %>%
    mutate(mutlist = ifelse(mutlist == "", "WT", mutlist)) 

dattoprint <- dattoprint %>% filter(!grepl("[A-Z][0-9]+X", mutlist))  %>% 
    mutate(mutlist = ifelse(!grepl("WT|M184", mutlist), "WT", mutlist))

dattoprint <- dattoprint %>% 
    group_by(loc, sampweek, mutlist) %>% summarize(freq = n()) %>% 
    arrange(loc, sampweek, desc(freq)) %>%
    spread(mutlist, freq, fill = 0)  %>%
    ungroup() %>% mutate(sampweek = 
              ifelse(sampweek == 12 | sampweek == 13, "5", 
              ifelse(sampweek == 15 | sampweek == 16, "20", 
              ifelse(sampweek == 20, "50", 
              ifelse(sampweek == 26, "100", NA )))))


dattoprint <- dattoprint %>% mutate(sampweek = paste0("samp", sampweek)) %>% 
    gather(mut, freq, -loc, -sampweek)   %>% 
    spread(sampweek, freq) %>% 
    select(loc, mut, `samp5`, `samp20`, `samp50`, `samp100`) %>%
    mutate(loc = gsub("RNA", "", loc)) %>% 
    arrange(loc, mut)


#########################################

#Exclude haplotypes that weren't seen before the third timepoint
exclude.haps <- paste((dattoprint %>% group_by(mut) %>% summarize(samp5 = sum(samp5), 
                                           samp20 = sum(samp20),
                                           samp50 = sum(samp50), 
                                           samp100 = sum(samp100)) %>% 
               filter(samp5 + samp20 == 0))$mut)
    
exc <- tbl_df(exclude.haps)
names(exc) <- c("key")

#This creates a map of the map-froms to the map-tos 
exc <- exc %>% mutate(simp = str_extract(gsub("-K249K", "", exc$key), "M184(I|V)(-N255N)?"))

#Replace all the late occurring haplotypes with their earlier occurring haplotypes
if(nrow(exc) >= 1){
dattoprint <- dattoprint %>% 
    mutate(mut = str_replace_all(mut, setNames(exc$simp, paste0("^", exc$key, "$")))) %>% 
    group_by(loc, mut) %>% 
    summarize(samp5 = sum(samp5), 
              samp20 = sum(samp20),
              samp50 = sum(samp50), 
              samp100 = sum(samp100))
}

#########################################

orderedtypes <- c("WT", (dattoprint %>% gather(sampweek, freq, -loc, -mut)  %>% 
    filter(mut != "WT") %>%
    group_by(mut) %>%
    summarize(n = sum(freq)) %>% arrange(desc(n)))$mut)


dattoprint <- dattoprint %>% mutate(mut = factor(mut, levels = orderedtypes)) %>% 
    group_by(loc) %>% arrange(loc, mut)

write.table(dattoprint, file = paste0("../output/shivdata-maf",maf,".txt"), quote= FALSE, row.names = FALSE, col.names = TRUE)

}
























