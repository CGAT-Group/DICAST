BEGIN{OFS=" "; print "id, transcript, strand, mate, start, end, range_needs_check"}
/@/ {
    gsub(/@/,""); 
    gsub("/","|");
    gsub(":","|");
    gsub("-","|");
	gsub(";","|")
    size = split($0,splitted,"|");
    if(size==8){
      print splitted[1], splitted[2], "*", "1", splitted[4], splitted[5], "FALSE"
      print splitted[1], splitted[2], "*", "2", splitted[7], splitted[8], "FALSE"
    }
    else{
      print splitted[1], splitted[2], "*", "1", "1", "76", "TRUE"
      print splitted[1], splitted[2], "*", "2", "1", "76", "TRUE"
    }
} 
    
