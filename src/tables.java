
public final class tables {

    static final public byte[] dx2 = new byte[128];
    static final public byte[] cdna = new byte[128];
    static final public byte[] cdnat2 = new byte[5];
    static final public byte[] cdn = new byte[128];
    static final public byte[] dx = new byte[128];
    
    static {
// "M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T and I"
//        a  CGT  c  AGT  a   c   g   ACT g   g   GT  t   AC  ATGC AG  GC  t  t   AGC AT   CT
//        a   b   c   d   e   f   g   h   i   j   k   l   m    n   r   s   t  u   v   w    y
//        97  98  99  100 101 102 103 104 105 106 107 108 109 110 114 115 116 117 118 119 121
        //antisense 
        cdnat2[0] = 1;
        cdnat2[1] = 0;
        cdnat2[2] = 3;
        cdnat2[3] = 2;
        cdnat2[4] = 4;

//' sp   -     a   b   c   d   e   f   g   h   i   j   k   l   m    n   r   s   t  u   v   w    y
//' 32  45    97  98  99  100 101 102 103 104 105 106 107 108 109 110 114 115 116 117 118 119 121
//' 0    0     1   2   3   4   5   6   7   8   9   10  11  12  13  14  18  19  20  21  22  23  25
        dx2[97] = 0;   //a =0
        dx2[101] = 0;  //a  
        dx2[108] = 1;  //t = 1            
        dx2[116] = 1;  //t 
        dx2[117] = 1;  //u  
        dx2[99] = 2;   //c = 2
        dx2[102] = 2;  //c  
        dx2[103] = 3;  //g = 3
        dx2[105] = 3;  //g 
        dx2[106] = 3;  //g    
        dx2[119] = 4;  //at  
        dx2[98] = 4;   //tgc   
        dx2[104] = 4;  //atc   
        dx2[121] = 4;  //tc        
        dx2[109] = 4;  //ac
        dx2[115] = 4;  //gc       
        dx2[107] = 4;  //gt   
        dx2[114] = 4;  //ag         
        dx2[118] = 4;  //agc   
        dx2[100] = 4;  //atg         
        dx2[110] = 4;  //atgc 

//M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T    
//a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   y
//97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 121            
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 99;   // I
        cdn[75] = 107;  // K
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y        
        cdn[97] = 97;   // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h
        cdn[105] = 99;   // i
        cdn[107] = 107;  // k
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y 

        cdna[65] = 84;  //A
        cdna[66] = 86;  //B
        cdna[67] = 71;  //C
        cdna[68] = 72;  //D
        cdna[71] = 67;  //G
        cdna[72] = 68;  //H
        cdna[73] = 71;  //I
        cdna[75] = 77;  //K
        cdna[77] = 75;  //M
        cdna[78] = 78;  //N
        cdna[82] = 89;  //R
        cdna[83] = 83;  //S
        cdna[84] = 65;  //T
        cdna[85] = 65;  //U
        cdna[86] = 66;  //V
        cdna[87] = 87;  //W
        cdna[89] = 82;  //Y
        
        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // i <- g
        cdna[107] = 109;// m <- k
        cdna[109] = 107;// k <- m
        cdna[110] = 110;// n <- n
        cdna[114] = 121;// y <- r
        cdna[115] = 115;// s <- s
        cdna[116] = 97; // a <- t
        cdna[117] = 97; // a <- u
        cdna[118] = 98; // b <- v
        cdna[119] = 119;// w <- w
        cdna[121] = 114;// r <- y
        
        //tb=dx
        dx[97] = 0;   //a =0
        dx[119] = 0;  //at                   
        dx[101] = 0;  //a  
        dx[108] = 1;  //t = 1            
        dx[116] = 1;  //t 
        dx[117] = 1;  //u  
        dx[98] = 2;   //tgc         
        dx[99] = 2;   //c = 2
        dx[102] = 2;  //c        
        dx[104] = 2;  //atc   
        dx[121] = 2;  //tc        
        dx[109] = 2;  //ac
        dx[115] = 2;  //gc        
        dx[103] = 3;  //g = 3
        dx[105] = 3;  //g 
        dx[106] = 3;  //g      
        dx[107] = 3;  //gt   
        dx[114] = 3;  //ag         
        dx[118] = 3;  //agc   
        dx[100] = 3;  //atg 
        dx[110] = 4;  //atgc 
    }
}
