	  x7  �   k820309    �          15.0        ��[                                                                                                           
       CABLE/cable_read.f90 READ_MODULE              REDISTR_I REDISTR_R REDISTR_RD REDISTR_R2 REDISTR_R2D gen@READPAR                                                    
                                                          
                            @                              
                            @                              
       LANDPT EXISTS LAND_X LAND_Y METGRID                                                        u #READPAR_I    #READPAR_R    #READPAR_RD    #READPAR_R2 *   #READPAR_R2D 6   #         @     @X                                             	   #READPAR_I%TRIM    #READPAR_I%PRESENT    #NCID    #PARNAME 	   #COMPLETESET 
   #VAR_I    #FILENAME    #NPATCH    #DIMSWITCH    #FROM_RESTART    #INPATCH                  @                                TRIM               @                                PRESENT           
  @                                                   
  @                              	                    1           
D                                 
                      
D @                                                                &                                                     
  @                                                  1           
                                                      
                                                     1           
 @                                                    
                                           #         @     @X                                             	   #READPAR_R%TRIM    #READPAR_R%PRESENT    #NCID    #PARNAME    #COMPLETESET    #VAR_R    #FILENAME    #NPATCH    #DIMSWITCH    #FROM_RESTART    #INPATCH                  @                                TRIM               @                                PRESENT           
  @                                                   
  @                                                  1           
D                                                       
D @                                                 	               &                                                     
  @                                                  1           
                                                      
                                                     1           
 @                                                    
                                           #         @     @X                                             	   #READPAR_RD%TRIM    #READPAR_RD%PRESENT    #READPAR_RD%REAL     #NCID !   #PARNAME "   #COMPLETESET #   #VAR_RD $   #FILENAME %   #NPATCH &   #DIMSWITCH '   #FROM_RESTART (   #INPATCH )                 @                                TRIM               @                                PRESENT               @            @                    REAL           
  @                              !                     
  @                              "                    1           
D                                 #                      
D                                $                   
               &                                                     
  @                              %                    1           
                                 &                     
                                 '                    1           
 @                               (                     
                                )           #         @     @X                             *                	   #READPAR_R2%TRIM +   #READPAR_R2%PRESENT ,   #NCID -   #PARNAME .   #COMPLETESET /   #VAR_R2 0   #FILENAME 1   #NPATCH 2   #DIMSWITCH 3   #FROM_RESTART 4   #INPATCH 5                 @                           +     TRIM               @                           ,     PRESENT           
  @                              -                     
  @                              .                    1           
D                                 /                      
D @                              0                   	               &                   &                                                     
  @                              1                    1           
                                 2                     
                                 3                    1           
 @                               4                     
                                5           #         @     @X                             6                	   #READPAR_R2D%TRIM 7   #READPAR_R2D%PRESENT 8   #READPAR_R2D%REAL 9   #NCID :   #PARNAME ;   #COMPLETESET <   #VAR_R2D =   #FILENAME >   #NPATCH ?   #DIMSWITCH @   #FROM_RESTART A   #INPATCH B                 @                           7     TRIM               @                           8     PRESENT               @            @              9     REAL           
  @                              :                     
  @                              ;                    1           
D                                 <                      
D                                =                   
               &                   &                                                     
  @                              >                    1           
                                 ?                     
                                 @                    1           
 @                               A                     
                                B                          @  @                          C     '                    #NAP D   #CSTART E   #CEND F   #ILAT G   #ILON H                �                              D                                �                              E                               �                              F                               �                              G                               �                              H                                 @  @                          I     '8                    #WIND J   #LWDOWN K   #CO2AIR L   #PSURF M   #SNOWF N   #AVPRECIP O   #LAI P   #LAI_T Q   #LAI_M R   #LAI_P S   #PARAMETERS T   #INITIAL U   #PATCH V   #LAIPATCH W                �                               J                                �                               K                               �                               L                               �                               M                               �                               N                               �                               O                               �                               P                               �                               Q                               �                               R             	                   �                               S     $       
                   �                               T     (                          �                               U     ,                          �                               V     0                          �                               W     4             %         @                                X                           #NCID Y   #NAME Z   #VARID [             
   @                              Y                     
   @                             Z                    1             @                              [                                                         \                                                       0%         @                                ]                          #NF90_INQUIRE_VARIABLE%PRESENT ^   #NF90_INQUIRE_VARIABLE%TRIM _   #NF90_INQUIRE_VARIABLE%SIZE `   #NCID a   #VARID b   #NAME c   #XTYPE d   #NDIMS e   #DIMIDS f   #NATTS g                 @                            ^     PRESENT               @                            _     TRIM               @                            `     SIZE           
   @                              a                     
   @                              b                       @                             c                     1             @                              d                        @                              e                       @                              f                                  &                                                       @                              g            #         @                                   h                    #INPATCH i   #NAP j   #IN_I k   #OUT_I l   #PARNAME n             
                                 i                    
                                 j                        p          5 � p        r i       5 � p        r i                              
                                 k                        p          5 � p        r i       5 � p        r i                              D                                l                         p          5 r m       5 r m                               
                                 n                    1 #         @                                   o                   #REDISTR_R%FLOAT p   #INPATCH q   #NAP r   #IN_R s   #OUT_R t   #PARNAME u                 @                           p     FLOAT           
                                 q                    
                                 r                        p          5 � p        r q       5 � p        r q                              
                                 s                    
    p          5 � p        r q       5 � p        r q                              D                                t                    
     p          5 r m       5 r m                               
                                 u                    1 #         @                                   v                   #REDISTR_RD%FLOAT w   #INPATCH x   #NAP y   #IN_RD z   #OUT_RD {   #PARNAME |                 @                           w     FLOAT           
                                 x                    
                                 y                         p          5 � p        r x       5 � p        r x                              
                                 z                    
 !   p          5 � p        r x       5 � p        r x                              D                                {                    
 "    p          5 r m       5 r m                               
                                 |                    1 #         @                                   }                   #REDISTR_R2%FLOAT ~   #INPATCH    #NAP �   #IN_R2 �   #OUT_R2 �   #PARNAME �   #DIM2 �                 @                           ~     FLOAT           
                                                     
                                 �                     #   p          5 � p        r        5 � p        r                               
                                 �                    
 $     p        5 � p        r    p          5 � p        r      5 � p        r �       5 � p        r      5 � p        r �                              D                                �                    
 %      p        5 r m   p          5 r m     5 � p        r �       5 r m     5 � p        r �                               
                                 �                    1           
                                 �           #         @                                   �                   #REDISTR_R2D%FLOAT �   #INPATCH �   #NAP �   #IN_R2D �   #OUT_R2D �   #PARNAME �   #DIM2 �                 @                           �     FLOAT           
                                 �                    
                                 �                     '   p          5 � p        r �       5 � p        r �                              
                                 �                    
 (     p        5 � p        r �   p          5 � p        r �     5 � p        r �       5 � p        r �     5 � p        r �                              D                                �                    
 )      p        5 r m   p          5 r m     5 � p        r �       5 r m     5 � p        r �                               
                                 �                    1           
                                 �                         @                             m               �   )      fn#fn !   �   R   b   uapp(READ_MODULE      @   J  ABORT_MODULE "   [  @   j  DEFINE_DIMENSIONS    �  @   J  NETCDF    �  d   J  IO_VARIABLES    ?  �       gen@READPAR    �  �      READPAR_I    �  =      READPAR_I%TRIM "   �  @      READPAR_I%PRESENT    9  @   a   READPAR_I%NCID "   y  L   a   READPAR_I%PARNAME &   �  @   a   READPAR_I%COMPLETESET       �   a   READPAR_I%VAR_I #   �  L   a   READPAR_I%FILENAME !   �  @   a   READPAR_I%NPATCH $     L   a   READPAR_I%DIMSWITCH '   i  @   a   READPAR_I%FROM_RESTART "   �  @   a   READPAR_I%INPATCH    �  �      READPAR_R    �  =      READPAR_R%TRIM "     @      READPAR_R%PRESENT    T  @   a   READPAR_R%NCID "   �  L   a   READPAR_R%PARNAME &   �  @   a   READPAR_R%COMPLETESET      	  �   a   READPAR_R%VAR_R #   �	  L   a   READPAR_R%FILENAME !   �	  @   a   READPAR_R%NPATCH $   8
  L   a   READPAR_R%DIMSWITCH '   �
  @   a   READPAR_R%FROM_RESTART "   �
  @   a   READPAR_R%INPATCH           READPAR_RD     
  =      READPAR_RD%TRIM #   G  @      READPAR_RD%PRESENT     �  =      READPAR_RD%REAL     �  @   a   READPAR_RD%NCID #     L   a   READPAR_RD%PARNAME '   P  @   a   READPAR_RD%COMPLETESET "   �  �   a   READPAR_RD%VAR_RD $     L   a   READPAR_RD%FILENAME "   h  @   a   READPAR_RD%NPATCH %   �  L   a   READPAR_RD%DIMSWITCH (   �  @   a   READPAR_RD%FROM_RESTART #   4  @   a   READPAR_RD%INPATCH    t  �      READPAR_R2     e  =      READPAR_R2%TRIM #   �  @      READPAR_R2%PRESENT     �  @   a   READPAR_R2%NCID #   "  L   a   READPAR_R2%PARNAME '   n  @   a   READPAR_R2%COMPLETESET "   �  �   a   READPAR_R2%VAR_R2 $   R  L   a   READPAR_R2%FILENAME "   �  @   a   READPAR_R2%NPATCH %   �  L   a   READPAR_R2%DIMSWITCH (   *  @   a   READPAR_R2%FROM_RESTART #   j  @   a   READPAR_R2%INPATCH    �  
     READPAR_R2D !   �  =      READPAR_R2D%TRIM $   �  @      READPAR_R2D%PRESENT !   1  =      READPAR_R2D%REAL !   n  @   a   READPAR_R2D%NCID $   �  L   a   READPAR_R2D%PARNAME (   �  @   a   READPAR_R2D%COMPLETESET $   :  �   a   READPAR_R2D%VAR_R2D %   �  L   a   READPAR_R2D%FILENAME #   *  @   a   READPAR_R2D%NPATCH &   j  L   a   READPAR_R2D%DIMSWITCH )   �  @   a   READPAR_R2D%FROM_RESTART $   �  @   a   READPAR_R2D%INPATCH '   6  �      LAND_TYPE+IO_VARIABLES +   �  H   a   LAND_TYPE%NAP+IO_VARIABLES .     H   a   LAND_TYPE%CSTART+IO_VARIABLES ,   I  H   a   LAND_TYPE%CEND+IO_VARIABLES ,   �  H   a   LAND_TYPE%ILAT+IO_VARIABLES ,   �  H   a   LAND_TYPE%ILON+IO_VARIABLES 0   !  �      INPUT_DETAILS_TYPE+IO_VARIABLES 5     H   a   INPUT_DETAILS_TYPE%WIND+IO_VARIABLES 7   _  H   a   INPUT_DETAILS_TYPE%LWDOWN+IO_VARIABLES 7   �  H   a   INPUT_DETAILS_TYPE%CO2AIR+IO_VARIABLES 6   �  H   a   INPUT_DETAILS_TYPE%PSURF+IO_VARIABLES 6   7  H   a   INPUT_DETAILS_TYPE%SNOWF+IO_VARIABLES 9     H   a   INPUT_DETAILS_TYPE%AVPRECIP+IO_VARIABLES 4   �  H   a   INPUT_DETAILS_TYPE%LAI+IO_VARIABLES 6     H   a   INPUT_DETAILS_TYPE%LAI_T+IO_VARIABLES 6   W  H   a   INPUT_DETAILS_TYPE%LAI_M+IO_VARIABLES 6   �  H   a   INPUT_DETAILS_TYPE%LAI_P+IO_VARIABLES ;   �  H   a   INPUT_DETAILS_TYPE%PARAMETERS+IO_VARIABLES 8   /  H   a   INPUT_DETAILS_TYPE%INITIAL+IO_VARIABLES 6   w  H   a   INPUT_DETAILS_TYPE%PATCH+IO_VARIABLES 9   �  H   a   INPUT_DETAILS_TYPE%LAIPATCH+IO_VARIABLES &     o       NF90_INQ_VARID+NETCDF +   v  @   e   NF90_INQ_VARID%NCID+NETCDF +   �  L   e   NF90_INQ_VARID%NAME+NETCDF ,      @   e   NF90_INQ_VARID%VARID+NETCDF "   B   q       NF90_NOERR+NETCDF -   �   �       NF90_INQUIRE_VARIABLE+NETCDF =   �!  @      NF90_INQUIRE_VARIABLE%PRESENT+NETCDF=PRESENT 7   �!  =      NF90_INQUIRE_VARIABLE%TRIM+NETCDF=TRIM 7   /"  =      NF90_INQUIRE_VARIABLE%SIZE+NETCDF=SIZE 2   l"  @   e   NF90_INQUIRE_VARIABLE%NCID+NETCDF 3   �"  @   e   NF90_INQUIRE_VARIABLE%VARID+NETCDF 2   �"  L   e   NF90_INQUIRE_VARIABLE%NAME+NETCDF 3   8#  @   e   NF90_INQUIRE_VARIABLE%XTYPE+NETCDF 3   x#  @   e   NF90_INQUIRE_VARIABLE%NDIMS+NETCDF 4   �#  �   e   NF90_INQUIRE_VARIABLE%DIMIDS+NETCDF 3   D$  @   e   NF90_INQUIRE_VARIABLE%NATTS+NETCDF    �$  �       REDISTR_I "   %  @   a   REDISTR_I%INPATCH    D%  �   a   REDISTR_I%NAP    �%  �   a   REDISTR_I%IN_I     �&  �   a   REDISTR_I%OUT_I "   @'  L   a   REDISTR_I%PARNAME    �'  �       REDISTR_R     !(  >      REDISTR_R%FLOAT "   _(  @   a   REDISTR_R%INPATCH    �(  �   a   REDISTR_R%NAP    S)  �   a   REDISTR_R%IN_R     *  �   a   REDISTR_R%OUT_R "   �*  L   a   REDISTR_R%PARNAME    �*  �       REDISTR_RD !   +  >      REDISTR_RD%FLOAT #   �+  @   a   REDISTR_RD%INPATCH    �+  �   a   REDISTR_RD%NAP !   �,  �   a   REDISTR_RD%IN_RD "   e-  �   a   REDISTR_RD%OUT_RD #   �-  L   a   REDISTR_RD%PARNAME    E.  �       REDISTR_R2 !   �.  >      REDISTR_R2%FLOAT #   %/  @   a   REDISTR_R2%INPATCH    e/  �   a   REDISTR_R2%NAP !   0  $  a   REDISTR_R2%IN_R2 "   =1  �   a   REDISTR_R2%OUT_R2 #   12  L   a   REDISTR_R2%PARNAME     }2  @   a   REDISTR_R2%DIM2    �2  �       REDISTR_R2D "   b3  >      REDISTR_R2D%FLOAT $   �3  @   a   REDISTR_R2D%INPATCH     �3  �   a   REDISTR_R2D%NAP #   �4  $  a   REDISTR_R2D%IN_R2D $   �5  �   a   REDISTR_R2D%OUT_R2D $   �6  L   a   REDISTR_R2D%PARNAME !   �6  @   a   REDISTR_R2D%DIM2 %   87  @      MP+DEFINE_DIMENSIONS 