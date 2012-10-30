# seq is the code need to be coded 'GATC' no 'N' 
# code.len is length of code, that's single nt code, di-nt code, tri-nt code ,etc
# seq is ordered as 'A' 'G' 'C' 'T'
# for example : 
# 	single nt : 
#		A : 000
#		C : 100
#		G : 010
#		T : 001
#	di-nt :
#		AA : 000,0000,0000,0000
#               AC : 100,0000,0000,0000
#               AG : 010,0000,0000,0000
#               AT : 001,0000,0000,0000
#               CA : 000,1000,0000,0000
#               CC : 000,0100,0000,0000
#               CG : 000,0010,0000,0000
#               CT : 000,0001,0000,0000
#               GA : 000,0000,1000,0000
#               GC : 000,0000,0100,0000
#               GG : 000,0000,0010,0000
#               GT : 000,0000,0001,0000
#               TA : 000,0000,0000,1000
#               TC : 000,0000,0000,0100
#               TG : 000,0000,0000,0010
#               TT : 000,0000,0000,0001
               
CodeSeq <- function(sequence,code.len=1)
{
	seq.len <- nchar(sequence[1])
	code <- matrix(0,length(sequence), (seq.len + 1 - code.len)*(4^code.len - 1) )
	idx.all <- matrix(0,length(sequence), (seq.len + 1 - code.len) )
	for (i in 1: (seq.len+1-code.len)){	
		code.cur <- matrix(0,length(sequence), 4^code.len - 1)	
		seq.cur <- substring(sequence,i,i+code.len-1)
		idx <- rep(0,length(sequence))
		for (j in 1:code.len){	
			idx[substring(seq.cur,j,j)=='A'] <- idx[substring(seq.cur,j,j)=='A'] + 0*4^(code.len-j)
			idx[substring(seq.cur,j,j)=='C'] <- idx[substring(seq.cur,j,j)=='C'] + 1*4^(code.len-j)
  			idx[substring(seq.cur,j,j)=='G'] <- idx[substring(seq.cur,j,j)=='G'] + 2*4^(code.len-j)
			idx[substring(seq.cur,j,j)=='T'] <- idx[substring(seq.cur,j,j)=='T'] + 3*4^(code.len-j)
		}	
		idx.all[,i] <- idx + (4^code.len-1)*(i-1)
		idx.all[idx==0,i] <- 0
	}
	
	for (i in 1:length(sequence)){
		code[i,idx.all[i,]]<- 1
		if (floor(i/100000)*100000 == i)
			print(i)
	}
	code
}



