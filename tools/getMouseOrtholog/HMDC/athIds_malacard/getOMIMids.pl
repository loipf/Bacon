#!/bin/perl

use IPC::System::Simple qw(system capture);

@disList=("ATH013","ART004","GNR003","CRB008","RNL008","ATH012","FGN001","ART021","ART140","CRN018","CRN300","END072","VSC007","HRT032","CRB009","ISC006","ISC004","AGN016","MYC007","HYP805","DBT009","CRB039","BDY004","BDY005","BDY019","BDY017","BDY015","BDY010","BDY011","BDY006","BDY012","HYP607","ART022","BLD140","CRT016","THR024","LPS004","RHM011","SYS001","HPT021","ANR040","KDN018","PRP080","CMP080","ANG054","PRD008","DBT026","CHR089","STR067","ATR011","URM002","ACT075","PRP027","CHL123","ALP042","HMZ003","OCL033","ATM095","HYP060","ALZ034","HYP595","ART016","FTT001","LPD008","ABD014","CRN030","CRT013","LVR013","HYP066","PNM007","CNN005","INT074","GLC003","HYP750","CNG034","PRT036","TNG002","PRD007","CLC006","ACQ007","SLP006","ABD017","ART106","BDY007","XNT003","RTN023","GLC008","AST005","SPN225","DMN002","VSC011","ENC018","STS003","URN009","LPD010","SPN051","INT007","FBR046","INF009","PSR001","NRP001","RHM027","ANT006","HYP724","VSC002","CLL015","CRN017","HYP555","MNT310","ATM052","KWS002","MCR115","PLY011","CNT035","HYP086","HYP266","INT002","LCT022","GT001","ERY003","ART138","HYP614","CRD132","HYP732","ATN004","RST023","DFF018","MCL042","OST002","BRC012","APN008","THL005","PSR002","BRG013","HYP190","TRN015","CYT008","MYP006","HYP014","OVR063","GRW007","LMB062","HPT025","BRN080","ENT004","ANG037","CTR119","MNS002","ATS285","HYP485","WRN001","ACR007","AMY004","ART023","HLC007","CRB011","HTC003","GLY013","JNT002","LYS012","END030","NNL002","HYP818","RTN022","CRN074","LPP002","RTR011","ACT088","BNM029","ART012","BNM022","DFC001","FML018","HPT001","ANR048","SKN016","PSD087","LNG099","VRL011","ABT001","RSP006","GST033","SLP005","HMR039","IMM136","MCR113","LPD015","HYP080","HYP739","ART115","HYP081","SKN027","PST011","PRD004","ART008","CRN019","MCR130","ARC023","MCR120","MCR133","ALR002","MDR006","LMN005","ANR038","CRT085","CRB087","BLD137","EPL164","MYC084","HPT016","SYS005","MXD005","ART067","ART101","GLL020","PNC044","VNW001","THR092","SCL052","NPH012","PRT013","HYP083","PTN001","THY032","CRD223","INF037","THY030","MYM013","HMC014","FML035","ART031","NPH010","LKD015","HYP734","HYP768","BNS003","CRT015","PDT035","EST007","HYP008","CRB088","CRN020","VSC018","TYM001","FML330","GPS001","ART110","FBR089","NST002","PSD059","BLD163","BLD138","BRS047","PRS040","CYS001","INF038","MYL005","PLM129","HPT073","ANX010","END044","CMM004","LYM118","SMT004","CHR012","HMP029","ALG028","TMP003","BRS051","LCL006","GLN010","INT068","MTH009","CRB175","PTT048","PRP029","HYP056","EXN002","MRB003","HMP007","GLM007","GST044","BRN019","PST028","HRD002","PRS047","PLM034","UND005","CRB172","MYC008","NNR004","THR015","ART017","THR100","IMM158","ANG015","PLM010","NPH009","RTN018","PYL005","ALP008","CRT049","ALK013","ALS001","THR004","PRM006","THR001","VNS003","FML012","OBS001","CHR005","CHL004","HYP806","CCN002","GNG011","FCT001","RNV001","DBT008","ACT078","ISC002","CHR008","BRN108","EPD070","THR035","RTN014","FSH001","HYP795","LYM024","IRN002","PRN037","SYS003","ADP001","HYP794","PLS029","HYP057","AND014","AMR003","LPD021","LPD035","ENC005","BLR024","CYT004","GND003","ACQ009","EHL081","ATY016","VRT003","ART010","NRG005","RNG031","LPD040","SLN001","STR035","LGP003","FLY003","ART153","ODN020","SBC025","ATR024","LPT011","RCT024","DMN012","OST012","MYL069","MLT020","GST053","INS024","PLM037","FLL037","PLY001","LKM002","MRF001","ALP046","MLD001","LVR012","CHG001","BRN028","NSP012","SCH015","WGN006","ESS003","JVN010","HRP006","DWN001","CNN003","ATN007","LYM017","HYD006","LYM007","ADL010","PLM164","CNJ013","HSH003","BRT054","SPN186","TYP007","TKY002","LPD012","THN009","CLC063","GST045","ADL030","SCR008","CRT072","RCT015","ALC007","MDD011","APP008","LPT001","FCT002","HYP020","CNT047","GNG013","SRC025","HPT003","KRT001","ACT210","MYS005","HMR004","CRD119","ALC004","INT066","SNS014","KRT019","ALC006","CMM005","CMR001","KBK002","FCT007","DFC004","MCR129","PRP030","NRV006","ACT074","SYN007","MCS002","TXC005","HMR003","MCP040","SPN052","CTR002","MNT002","IDP011","3HY005","ALL010","ACN002","PLY018","DBT010","INF032","LGG001","ALZ049","PRC002","BCT007","PLS007","CHL067","HYP069","INT030","PLY019","FLR002","CMP010","OVR049","PLM012","END033","ACT105","SJG008","GLY060","LKD001","CYS018","CNS004","ERD001","SLC006","EXT034","BRN106","ISL003","PTT009","INF034","BCT002","DBN001","BNF002","PRN019","LYM027","RTN016","IGR001","PLS006","PLV003","MGR003","NPH018","EXF001","CYS010","MCN007","IRN001","BRN022","AGG001","AMN001","SPS003","PRC013","SYP003","TST014","TRC086","RHM028","SCK005","IMP005","LKN025","CLC001","SCK002","HYP005","HPT009","ANT034","VRL007","CTN003","TRN034","DDN001","FBR032","STM007","PST095","OVR046","HYP068","VTM028","ILT001","PLR008","FNC043","PRC012","DBT004","PLY017","NPH091","FML089","GSG001","ATR005","END046","OPN001","INT067","PST041","STR072","HYP043","PPL049","DNT012","ACT017","HMC002","CLD007","PRM108","SCL015","HMS001","ACT068","URT031","IDP024","CRV039","VRL012","PRD040","MTB004","ACT164","DBT006","INH020","PTT037","PLM017","HYP740","VND001","PRP036","AMY009","DWR001","CNT046","DST006","CHR431","HRM006","DNG001","NM001","HYP458","NTR018","ACD008","SPN119","HYP085","FTL021","NSP002","MTR002","DDN006","RTN021","PTN014","CRD137","OBS037","HRN001","CNT015","ACQ022","SXL003","ANL018","MSN001","DFF003","MYX004","ALP061","CRS001","MDS022","CRT009","KRT008","TTR005","MNN034","OST062","NRN002","CHR579","IDP033","DFF006","ANC002","HYP781","CLP006","MYL013","SPS219","ART018","CRB086","ERL029","TND005","BCL014","HMF004","LTH045","CHR413","SHW001","RHM014","ALX002","PTC002","ATM022","FNS001","CNT028","PST055","STC004","CHR096","LNG057","MRN002","ACL001","PLY112","ALN001","IMM001","MLT002","LKS001","SKN006","AMN002","GRN009","HYP064","RNL097","LSB001","ALC010","BLD051","XNT009","CRT004","ALP041","PRG003","HYP396","SLT014","SBC035","CMP028","CRB194","CRC039","AML044","HNT004","LPT014","CMP090","HPT020","VSC008","DNS007","APL017","MLT151","INF013","PHG002","HRP025","HYP072","ACT036","DYS013","HYP114","EPT021","HNT011","SPL006","LYM045","MNG003","SPN331","LKC003","HMP001","BCT020","TRC078","SCR039","ACT209","EXT035","MLT009","AGR019","ANG066","CTR102","HRT008","CRN006","IDP041","RNL001","SBD001","BRK012","SBN001","CHR013","PRL001","ATM068","PYS001","PLY105","LPG001","PLT031","ANR010","TRT003","MLG042","PRM092","CRB031","SCL022","ART028","SPR024","CLS047","BRN041","FTT008","PLM064","INT053","XLN206","HYP002","RDN004","SKN010","MLG003","FML206","GRW032","ART134","EPS033","MRB006","PGT011","LPD028","3P2001","HYP163","ART071","LTH004","GST056","INB001","VSC001","ATY022","AMN012","CNG336","ANR037","CRT043","PLS026","VSC009","PRP050");
#@disList=("ATH013","ART004","GNR003");

for my $n (@disList) {
	# incapsula crawler
	my $url = capture("python3", "malacard_incapsula.py", $n);
	#my $url = `curl -s "https://www.malacards.org/card/$n"`; # incapsula protected
	print $n.",";
	
	# get OMIM id
	if( $url =~ /omim.org\/entry\/\d+"\starget="_blank">(\d{4,})<\/a>\s+<\/div>/) {
		print $1.",";
	} else {
		print ",";
	}

	# get DO id
	if( $url =~ /disease-ontology.org\/term\/DOID:\d+"\starget="_blank">(DOID:\d{4,})<\/a>\s+<\/div>/) {
		print $1;
	} 
	print "\n"
}










