package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	. "github.com/meoow/gomiscutils"
	"os"
	"regexp"
	"strconv"
	"strings"
)

const (
	VERSION string = "0.3.0a"
	GRCh    string = "GRCh37"
	CHR7    string = "CRA_TCAGchr7v2"
)

type (
	SEQGENEDB [24]CONTIG         // { 0, 1, 2, ... 21 , X, Y}
	CONTIG    map[string]*ORIENT //map[contig_id](map[gene_id]*GENEINFO)
	ORIENT    struct {
		PLUS  GENES
		MINUS GENES
	}
	GENES    map[uint32]*GENEINFO //map[gene_ID]*GENEINFO
	GENEINFO struct {
		GENENAME string
		P5       uint32
		P3       uint32
		MAP      string
	}
	RS2GENE_DIST_INFO_LIST struct {
		DISTINFO []*GENE_DIST_INFO
	}

	GENE_DIST_INFO struct {
		GENEID   uint32
		GENENAME string
		DIST     uint32
	}
)

var (
	chr7genes             = make(map[uint32]struct{}, 2000)
	groupLabelPattern     = regexp.MustCompile(`^GRCh37\.p10`)
	chrNumPattern         = regexp.MustCompile(`^([1-9][0-9]?|[Xx]|[Yy])(\|.+)?$`)
	chr7groupLabelPattern = regexp.MustCompile(`^CRA_TCAGchr7v2`)
	OrientationSignal     = map[int]string{0: "+", 1: "-"}
	dbgenefile            string
	snplist               string
	excludedsnplist       string
	version               bool
	extractGene           bool
	dist_th               string
	chrDist               bool
	Dist                  int64
)

const (
	field_isGene     = 11
	field_chrNum     = 1
	field_groupLabel = 12
	field_contig     = 5
	field_orient     = 8
	field_geneID     = 10
	field_geneName   = 9
	field_ctgStart   = 6
	field_ctgEnd     = 7
	field_chrStart   = 2
	field_chrEnd     = 3
)

var (
	field_cStart int
	field_cEnd   int
)

func init() {
	flag.StringVar(&dbgenefile, "db", "GENE.GZ", "Database file from dbSNP (gz compressed)")
	flag.StringVar(&snplist, "rs", "RSLIST.TXT", "One rs info per line listed file.\n        RS  CHR  CONTIG  GENE_POS  CHROME_POS(no use for now)")
	flag.BoolVar(&chrDist, "c", false, "Use CHR dist instead of GENE dist")
	flag.BoolVar(&extractGene, "p", false, "Just print gene db")
	flag.StringVar(&excludedsnplist, "e", "", "A file with excluded SNP list")
	flag.StringVar(&dist_th, "d", "-1", "Use a constant disance as threashold, may use (k)")
	flag.BoolVar(&version, "v", false, "Show version")
}

func main() {
	flag.Parse()
	if version {
		fmt.Println(VERSION)
		os.Exit(1)
	}
	//	if flag.NFlag() <= 2 {
	//		flag.PrintDefaults()
	//		os.Exit(0)
	//	}
	dist_th = strings.ToLower(dist_th)
	var multiply = 1
	if strings.HasSuffix(dist_th, "k") {
		dist_th = strings.TrimSuffix(dist_th, "k")
		multiply = 1000
	}
	dist, err := strconv.ParseInt(dist_th, 10, 32)
	Die(err)
	dist *= int64(multiply)
	Dist = dist
	//os.Stderr.WriteString(fmt.Sprintf("%d", Dist))

	if chrDist {
		field_cStart = field_chrStart
		field_cEnd = field_chrEnd
	} else {
		field_cStart = field_ctgStart
		field_cEnd = field_ctgEnd
	}

	var seqgendb SEQGENEDB
	for i := range seqgendb {
		seqgendb[i] = make(CONTIG, 20)
	}
	(&seqgendb).ReadDBFile(dbgenefile)
	(&seqgendb).ReadDBFileForChr7(dbgenefile)
	if extractGene {
		(&seqgendb).PrintDB()
		os.Exit(0)
	}
	(&seqgendb).ParseRSList(snplist, excludedsnplist)
}

func (db *SEQGENEDB) ReadDBFile(pathOfFile string) {
	var (
		f                *os.File
		err              error
		fReader          *bufio.Reader
		gzReader         *gzip.Reader
		contig           string
		chrnum           uint32
		genename         string
		geneid           uint32
		ctgOrient        string
		ctgPos1, ctgPos2 uint32
		//orient uint8
	)

	f, err = os.Open(pathOfFile)
	Die(err)
	defer f.Close()

	if isGzip() {
		gzReader, err = gzip.NewReader(f)
		Die(err)
		fReader = bufio.NewReader(gzReader)
	} else {
		fReader = bufio.NewReader(f)
	}
	lines := Readline(fReader)
	for line := range lines {
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")

		//		if len(fields[14]) < 4 || fields[14][0:4] != "best" {
		//			continue }
		if fields[field_isGene] != "GENE" {
			continue
		}
		if !chrNumPattern.MatchString(fields[field_chrNum]) {
			continue
		}
		if !(groupLabelPattern.MatchString(fields[field_groupLabel])) {
			continue
		}

		contig, ctgOrient, genename, chrnum, ctgPos1, ctgPos2, geneid = parseDBFieldsInfo(fields)

		if _, ok := db[chrnum][contig]; !ok {
			db[chrnum][contig] = &ORIENT{make(GENES, 30), make(GENES, 30)}
		}
		switch ctgOrient {
		case "+":
			if _, ok := db[chrnum][contig].PLUS[geneid]; !ok {
				if chrnum+1 == 7 {
					chr7genes[geneid] = struct{}{}
				}
				db[chrnum][contig].PLUS[geneid] = &GENEINFO{genename, 0, 0, GRCh}
			}
			db[chrnum][contig].PLUS[geneid].P5 = ctgPos1
			db[chrnum][contig].PLUS[geneid].P3 = ctgPos2
		case "-":
			if _, ok := db[chrnum][contig].MINUS[geneid]; !ok {
				if chrnum+1 == 7 {
					chr7genes[geneid] = struct{}{}
				}
				db[chrnum][contig].MINUS[geneid] = &GENEINFO{genename, 0, 0, GRCh}
			}
			db[chrnum][contig].MINUS[geneid].P3 = ctgPos1
			db[chrnum][contig].MINUS[geneid].P5 = ctgPos2
		default:
			panic("ORIENT TYPE ERROR: unknown orientation signal.")
		}
	}
}

func (db *SEQGENEDB) ReadDBFileForChr7(pathOfFile string) {
	var (
		f                *os.File
		err              error
		fReader          *bufio.Reader
		gzReader         *gzip.Reader
		contig           string
		chrnum           uint32
		genename         string
		geneid           uint32
		ctgOrient        string
		ctgPos1, ctgPos2 uint32
		//var orient uint8
	)

	f, err = os.Open(pathOfFile)
	Die(err)
	defer f.Close()

	if isGzip() {
		gzReader, err = gzip.NewReader(f)
		Die(err)
		fReader = bufio.NewReader(gzReader)
	} else {
		fReader = bufio.NewReader(f)
	}
	lines := Readline(fReader)
	for line := range lines {
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")

		if !chrNumPattern.MatchString(fields[field_chrNum]) {
			continue
		}

		if parseChrNum(fields[field_chrNum]) != "7" {
			continue
		}
		//		if len(fields[14]) < 4 || fields[14][0:4] != "best" {
		//			continue }
		if fields[field_isGene] != "GENE" {
			continue
		}
		if !chrNumPattern.MatchString(fields[field_chrNum]) {
			continue
		}
		if !(chr7groupLabelPattern.MatchString(fields[field_groupLabel])) {
			continue
		}

		contig, ctgOrient, genename, chrnum, ctgPos1, ctgPos2, geneid = parseDBFieldsInfo(fields)
		if _, ok := chr7genes[geneid]; ok {
			continue
		}
		if _, ok := db[chrnum][contig]; !ok {
			db[chrnum][contig] = &ORIENT{make(GENES, 30), make(GENES, 30)}
		}
		switch ctgOrient {
		case "+":
			if _, ok := db[chrnum][contig].PLUS[geneid]; !ok {
				//os.Stderr.WriteString("chr7 gene mapped\n")
				db[chrnum][contig].PLUS[geneid] = &GENEINFO{genename, 0, 0, CHR7}
				db[chrnum][contig].PLUS[geneid].P5 = ctgPos1
				db[chrnum][contig].PLUS[geneid].P3 = ctgPos2
			}
		case "-":
			if _, ok := db[chrnum][contig].MINUS[geneid]; !ok {
				//os.Stderr.WriteString("chr7 gene mapped\n")
				db[chrnum][contig].MINUS[geneid] = &GENEINFO{genename, 0, 0, CHR7}
				db[chrnum][contig].MINUS[geneid].P3 = ctgPos1
				db[chrnum][contig].MINUS[geneid].P5 = ctgPos2
			}
		default:
			panic("ORIENT TYPE ERROR: unknown orientation signal.")
		}
	}
}

func parseDBFieldsInfo(f []string) (contig, orient, gname string, chrnum, ctg1, ctg2, gid uint32) {
	var chrnum_string string
	switch parseChrNum(f[field_chrNum]) {
	case "x":
		fallthrough
	case "X":
		chrnum_string = "23"
	case "y":
		fallthrough
	case "Y":
		chrnum_string = "24"
	default:
		chrnum_string = parseChrNum(f[field_chrNum])
	}
	contig = f[field_contig]
	orient = f[field_orient]
	chrnum = uint32(MustParseUint(chrnum_string, 10, 32)) - 1
	ctg1 = uint32(MustParseUint(f[field_cStart], 10, 32))
	ctg2 = uint32(MustParseUint(f[field_cEnd], 10, 32))
	gname = f[field_geneName]
	gid = uint32(MustParseUint(f[field_geneID][7:], 10, 32))
	return
}

func parseChrNum(num string) string {
	return chrNumPattern.FindStringSubmatch(num)[1]
}

func (db *SEQGENEDB) PrintDB() {
	for i, contigs := range db {
		for contig, orient := range contigs {
			for geneid, geneinfo := range orient.PLUS {
				fmt.Printf("%2d%16s%12d%18s%10d%10d%3s%16s\n", i+1, contig, geneid, geneinfo.GENENAME, geneinfo.P5, geneinfo.P3, "+", geneinfo.MAP)
			}
			for geneid, geneinfo := range orient.MINUS {
				fmt.Printf("%2d%16s%12d%18s%10d%10d%3s%16s\n", i+1, contig, geneid, geneinfo.GENENAME, geneinfo.P3, geneinfo.P5, "-", geneinfo.MAP)
			}
		}
	}
}

func isGzip() bool {
	return true
}

func (db *SEQGENEDB) ParseRSList(pathOfFile string, exclFile string) {
	var f *os.File
	var err error
	var fReader *bufio.Reader
	var rsnum uint32
	var contig string
	var chrnum uint32
	var pos uint32
	var p3dist, p5dist uint32
	var lostnfound = make([]uint32, 0, 1000)
	var exclSNP = parseExcluRSList(exclFile)
	rs2genP5DistList := &RS2GENE_DIST_INFO_LIST{make([]*GENE_DIST_INFO, 0, 2)}
	rs2genP3DistList := &RS2GENE_DIST_INFO_LIST{make([]*GENE_DIST_INFO, 0, 2)}

	printOutR2GMap := func(rsnum uint32, r2gdl *RS2GENE_DIST_INFO_LIST, p uint8) {
		for {
			i := r2gdl.Pop()
			if i != nil {
				fmt.Printf("%d\t%d\t%s\t%d\t%d\n", rsnum, i.GENEID,
					i.GENENAME, i.DIST, p)
			} else {
				break
			}
		}
	}
	printLostAndFound := func(list []uint32) {
		os.Stderr.WriteString("# The Following SNP ain't found in database:\n")
		for _, v := range list {
			os.Stderr.WriteString(fmt.Sprintf("%d\n", v))
		}
	}
	reachBottom := func(a, b uint32) uint32 {
		if a > b {
			return a - b
		} else {
			return 0
		}
		panic("")
	}
	f, err = os.Open(pathOfFile)
	Die(err)
	defer f.Close()

	fReader = bufio.NewReader(f)
	lines := Readline(fReader)
	for line := range lines {
		fields := strings.Split(TrimNewLine(line), "\t")
		switch fields[1] {
		case "x":
			fallthrough
		case "X":
			fields[1] = "23"
		case "y":
			fallthrough
		case "Y":
			fields[1] = "24"
		}
		if !chrNumPattern.MatchString(fields[1]) {
			continue
		}
		rsnum = uint32(MustParseUint(fields[0], 10, 32))
		if _, ok := exclSNP[rsnum]; ok {
			continue
		}
		chrnum = uint32(MustParseUint(fields[1], 10, 32)) - 1
		contig = fields[2]
		pos = uint32(MustParseUint(fields[3], 10, 32))
		rs2genP5DistList.RemoveAll()
		rs2genP3DistList.RemoveAll()
		//		rs2genP5DistList.Push(&GENE_DIST_INFO{0, "", math.MaxUint32})
		//		rs2genP3DistList.Push(&GENE_DIST_INFO{0, "", math.MaxUint32})
		if _, ok := db[chrnum][contig]; !ok {
			lostnfound = append(lostnfound, rsnum)
			continue
		}
		for gid, ginfo := range db[chrnum][contig].PLUS {
			if isBetweenRange(reachBottom(ginfo.P5, 2000), ginfo.P3+500, pos) {
				p5dist = 0
				p3dist = 0
			} else {
				p5dist = simpleSubstractAbs(ginfo.P5, pos)
				p3dist = simpleSubstractAbs(ginfo.P3, pos)
			}
			rs2genP5DistList.Push(&GENE_DIST_INFO{gid, ginfo.GENENAME, p5dist})
			rs2genP3DistList.Push(&GENE_DIST_INFO{gid, ginfo.GENENAME, p3dist})
		}
		for gid, ginfo := range db[chrnum][contig].MINUS {
			if isBetweenRange(reachBottom(ginfo.P3, 500), ginfo.P5+2000, pos) {
				p5dist = 0
				p3dist = 0
			} else {
				p5dist = simpleSubstractAbs(ginfo.P5, pos)
				p3dist = simpleSubstractAbs(ginfo.P3, pos)
			}
			rs2genP5DistList.Push(&GENE_DIST_INFO{gid, ginfo.GENENAME, p5dist})
			rs2genP3DistList.Push(&GENE_DIST_INFO{gid, ginfo.GENENAME, p3dist})
		}
		printOutR2GMap(rsnum, rs2genP5DistList, 5)
		printOutR2GMap(rsnum, rs2genP3DistList, 3)
	}
	printLostAndFound(lostnfound)
}

func (info *RS2GENE_DIST_INFO_LIST) Push(r *GENE_DIST_INFO) {
	if Dist < 0 {
		switch {
		case len(info.DISTINFO) < 0:
			goto APPEND
		case info.DISTINFO[0].DIST > r.DIST:
			info.RemoveAll()
			fallthrough
		case info.DISTINFO[0].DIST == r.DIST:
			goto APPEND
		default:
			return
		}
	} else if int64(r.DIST) > Dist {
		return
	}
APPEND:
	info.DISTINFO = append(info.DISTINFO, r)
	return
}

func (info *RS2GENE_DIST_INFO_LIST) RemoveAll() {
	info.DISTINFO = info.DISTINFO[:0]
}

func (info *RS2GENE_DIST_INFO_LIST) String() string {
	return fmt.Sprintf("%v", info.DISTINFO)
}

func (info *RS2GENE_DIST_INFO_LIST) Dist() uint32 {
	return info.DISTINFO[0].DIST
}

func (info *RS2GENE_DIST_INFO_LIST) InGeneID(geneid uint32) bool {
	for _, i := range info.DISTINFO {
		if geneid == i.GENEID {
			return true
		}
	}
	return false
}

func (info *RS2GENE_DIST_INFO_LIST) Pop() *GENE_DIST_INFO {
	l := len(info.DISTINFO)
	if l == 0 {
		return nil
	} else {
		last := info.DISTINFO[l-1]
		info.DISTINFO = info.DISTINFO[:l-1]
		return last
	}
}

func simpleSubstractAbs(a, b uint32) uint32 {
	switch {
	case a >= b:
		return a - b
	default:
		return b - a
	}
	panic("")
}

func isBetweenRange(s, e, i uint32) bool {
	if s > e {
		s, e = e, s
	}
	if s <= i && e >= i {
		return true
	} else {
		return false
	}
}

func parseExcluRSList(pathOfFile string) map[uint32]struct{} {
	exclSNP := make(map[uint32]struct{})
	if pathOfFile != "" {
		f, e := os.Open(pathOfFile)
		Die(e)
		defer f.Close()
		fr := bufio.NewReader(f)
		lines := Readline(fr)
		for line := range lines {
			line = TrimNewLine(line)
			if line[0:2] == "rs" {
				line = line[2:]
			}
			exclSNP[uint32(MustParseUint(line, 10, 32))] = struct{}{}
		}
	}
	return exclSNP
}
