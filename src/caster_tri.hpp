#ifndef CASTER_TRI_HPP
#define CASTER_TRI_HPP

#include "stepwise_colorable.hpp"
#include "alignment_utilities.hpp"

namespace caster_tri{

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;
using std::unordered_map;

// +++ CASTER_TRI: static cache of priority taxa for taxonOrderPrioritizing (set during read()) +++
inline std::vector<size_t> s_priorityTaxa;
// +++ END CASTER_TRI +++
	
namespace DriverHelper {
	template<typename DataClasses> DataClasses read();
};

template<class Attributes> concept STEPWISE_COLOR_ATTRIBUTES = requires
{
    requires std::integral<typename Attributes::score_t> || std::floating_point<typename Attributes::score_t>;
	requires std::integral<typename Attributes::cnt_t> || std::floating_point<typename Attributes::cnt_t>;
	requires std::integral<typename Attributes::cnt4_t> || std::floating_point<typename Attributes::cnt4_t>;
	requires std::integral<typename Attributes::index_t>;
	{ Attributes::ZERO } -> std::convertible_to<typename Attributes::score_t>;
	{ Attributes::EPSILON } -> std::convertible_to<typename Attributes::score_t>;
};

template<typename cnt_taxon_type = unsigned char, typename cnt_type = unsigned short> struct StepwiseColorDefaultAttributes {
	using score_t = double;
	using cnt_taxon_t = cnt_taxon_type;
	using cnt_t = cnt_type;
	using cnt4_t = std::conditional_t<std::same_as<cnt_type, unsigned char>, unsigned int, unsigned long long>;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = 1e-3;
};

ChangeLog logColor("Color",
	"2026-02-02", "Chao Zhang", "Supporting quadripartiton", "minor",
	// +++ CASTER_TRI +++
	"2026-06-12", "Zuizhi Chen", "CASTER_TRI: only score quartets containing reference species", "minor",
	"2026-06-15", "Zuizhi Chen", "CASTER_TRI: adapt to TAXON_ORDER_PRIORITIZING interface", "minor");
	// +++ END CASTER_TRI +++

template<STEPWISE_COLOR_ATTRIBUTES Attributes> class Color{
	using cnt_taxon_t = Attributes::cnt_taxon_t;
	using cnt_t = Attributes::cnt_t;
	using cnt4_t = Attributes::cnt4_t;
	using index_t = Attributes::index_t;

public:
	using score_t = Attributes::score_t;
	static inline bool constexpr IS_ROOTED = false;
	static inline score_t constexpr ZERO = Attributes::ZERO;
	static inline score_t constexpr EPSILON = Attributes::EPSILON;
	
	struct SharedConstData{
		using ParentClass = Color<Attributes>;

		struct Element{
			index_t iGenomePosBegin = 0;
			index_t nPos = 0;
			vector<vector<array<cnt_taxon_t, 4> > > cnts; // cnts[iRow][iPos][iNucleotide] -> count
			vector<index_t> taxon2row; // taxon2row[iTaxon] -> iRow in cnts
			array<score_t, 4> eqFreqs{}; // eqFreqs[iNucleotide]
			// +++ CASTER_TRI: which reference species this element was aligned to +++
			size_t iReferenceTaxonId = (size_t)-1;
			// +++ END CASTER_TRI +++

			bool hasTaxon(size_t iTaxon) const noexcept{
				return iTaxon < taxon2row.size() && taxon2row[iTaxon] != -1;
			}
		};

		vector<Element> elements;
		index_t nGenomePos = 0;
		// +++ CASTER_TRI: all reference species taxon IDs used for priority ordering +++
		vector<size_t> priorityTaxa;
		// +++ END CASTER_TRI +++

		size_t nElements() const noexcept { return elements.size(); }
    };

	// +++ CASTER_TRI: static method for TAXON_ORDER_PRIORITIZING concept +++
	// Moves all reference species to the front of taxonOrder (preserving their shuffled relative order)
	static void taxonOrderPrioritizing(std::vector<size_t>& taxonOrder) noexcept {
		if (s_priorityTaxa.empty()) return;
		std::vector<size_t> priority, others;
		for (size_t t : taxonOrder) {
			if (std::find(s_priorityTaxa.begin(), s_priorityTaxa.end(), t) != s_priorityTaxa.end())
				priority.push_back(t);
			else
				others.push_back(t);
		}
		size_t i = 0;
		for (size_t t : priority) taxonOrder[i++] = t;
		for (size_t t : others) taxonOrder[i++] = t;
	}
	// +++ END CASTER_TRI +++

private:
	SharedConstData const& sharedConstData;
    vector<array<array<cnt_t, 4>, 4> > colorCnts; // colorCnts[iGenomePos][iColor][iNucleotide] -> count
	// +++ CASTER_TRI: reference species counts for R-subtraction; mutually exclusive w.r.t element refs +++
	vector<array<array<cnt_t, 4>, 4> > rColorCnts; // rColorCnts[iGenomePos][iColor][iNucleotide] -> count of ref only
	// +++ END CASTER_TRI +++

	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];
		if (!element.hasTaxon(iTaxon)) return;
		index_t iRow = element.taxon2row[iTaxon];
		index_t iPosBegin = element.iGenomePosBegin;
		for (index_t iPos : iota((index_t)0, element.nPos)){
			for (index_t iNucleotide : iota((index_t)0, (index_t)4)) {
				cnt_t& colorCnt = colorCnts[iPosBegin + iPos][iColor][iNucleotide];
				cnt_t cnt = element.cnts[iRow][iPos][iNucleotide];
				if constexpr (isSet) colorCnt += cnt;
				else colorCnt -= cnt;
			}
		}
		// +++ CASTER_TRI: if this taxon is the ref species for this element, also track its counts +++
		if (iTaxon == element.iReferenceTaxonId) {
			for (index_t iPos : iota((index_t)0, element.nPos)){
				for (index_t iNucleotide : iota((index_t)0, (index_t)4)) {
					cnt_t& rCnt = rColorCnts[iPosBegin + iPos][iColor][iNucleotide];
					cnt_t cnt = element.cnts[iRow][iPos][iNucleotide];
					if constexpr (isSet) rCnt += cnt;
					else rCnt -= cnt;
				}
			}
		}
		// +++ END CASTER_TRI +++
	}
	
	inline static cnt4_t XXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t y0, cnt4_t y1, cnt4_t y2) noexcept{
		return x0 * (x0 - 1) * y1 * y2 + x1 * (x1 - 1) * y2 * y0 + x2 * (x2 - 1) * y0 * y1
			 + y0 * (y0 - 1) * x1 * x2 + y1 * (y1 - 1) * x2 * x0 + y2 * (y2 - 1) * x0 * x1;
	}

	// +++ CASTER_TRI: internal score helper operating on cnt4_t values +++
	inline static score_t scorePosFromValues(
		cnt4_t a0, cnt4_t c0, cnt4_t g0, cnt4_t t0,
		cnt4_t a1, cnt4_t c1, cnt4_t g1, cnt4_t t1,
		cnt4_t a2, cnt4_t c2, cnt4_t g2, cnt4_t t2,
		array<score_t, 4> const &pi) noexcept{

		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt4_t const r0 = a0 + g0, y0 = c0 + t0;
		cnt4_t const r1 = a1 + g1, y1 = c1 + t1;
		cnt4_t const r2 = a2 + g2, y2 = c2 + t2;

		cnt4_t const rryy = XXYY(r0, r1, r2, y0, y1, y2);

		cnt4_t const aayy = XXYY(a0, a1, a2, y0, y1, y2);
		cnt4_t const ggyy = XXYY(g0, g1, g2, y0, y1, y2);
		cnt4_t const rrcc = XXYY(r0, r1, r2, c0, c1, c2);
		cnt4_t const rrtt = XXYY(r0, r1, r2, t0, t1, t2);

		cnt4_t const aacc = XXYY(a0, a1, a2, c0, c1, c2);
		cnt4_t const aatt = XXYY(a0, a1, a2, t0, t1, t2);
		cnt4_t const ggcc = XXYY(g0, g1, g2, c0, c1, c2);
		cnt4_t const ggtt = XXYY(g0, g1, g2, t0, t1, t2);

		return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
			 + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
	}

	inline static score_t scorePos(array<array<cnt_t, 4>, 4> const &cnt, array<score_t, 4> const &pi) noexcept{
		return scorePosFromValues(
			(cnt4_t)cnt[0][0], (cnt4_t)cnt[0][1], (cnt4_t)cnt[0][2], (cnt4_t)cnt[0][3],
			(cnt4_t)cnt[1][0], (cnt4_t)cnt[1][1], (cnt4_t)cnt[1][2], (cnt4_t)cnt[1][3],
			(cnt4_t)cnt[2][0], (cnt4_t)cnt[2][1], (cnt4_t)cnt[2][2], (cnt4_t)cnt[2][3],
			pi);
	}

	// +++ CASTER_TRI: scorePos with R-subtraction (only count quartets containing R) +++
	inline static score_t scorePosMinusR(array<array<cnt_t, 4>, 4> const &cnt, array<array<cnt_t, 4>, 4> const &rCnt, array<score_t, 4> const &pi) noexcept{
		score_t all = scorePosFromValues(
			(cnt4_t)cnt[0][0], (cnt4_t)cnt[0][1], (cnt4_t)cnt[0][2], (cnt4_t)cnt[0][3],
			(cnt4_t)cnt[1][0], (cnt4_t)cnt[1][1], (cnt4_t)cnt[1][2], (cnt4_t)cnt[1][3],
			(cnt4_t)cnt[2][0], (cnt4_t)cnt[2][1], (cnt4_t)cnt[2][2], (cnt4_t)cnt[2][3],
			pi);
		score_t minusR = scorePosFromValues(
			(cnt4_t)cnt[0][0] - (cnt4_t)rCnt[0][0], (cnt4_t)cnt[0][1] - (cnt4_t)rCnt[0][1], (cnt4_t)cnt[0][2] - (cnt4_t)rCnt[0][2], (cnt4_t)cnt[0][3] - (cnt4_t)rCnt[0][3],
			(cnt4_t)cnt[1][0] - (cnt4_t)rCnt[1][0], (cnt4_t)cnt[1][1] - (cnt4_t)rCnt[1][1], (cnt4_t)cnt[1][2] - (cnt4_t)rCnt[1][2], (cnt4_t)cnt[1][3] - (cnt4_t)rCnt[1][3],
			(cnt4_t)cnt[2][0] - (cnt4_t)rCnt[2][0], (cnt4_t)cnt[2][1] - (cnt4_t)rCnt[2][1], (cnt4_t)cnt[2][2] - (cnt4_t)rCnt[2][2], (cnt4_t)cnt[2][3] - (cnt4_t)rCnt[2][3],
			pi);
		return all - minusR;
	}
	// +++ END CASTER_TRI +++
	
	inline static cnt4_t quadXXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t x3, cnt4_t y0, cnt4_t y1, cnt4_t y2, cnt4_t y3) noexcept {
		return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
	}

	static score_t quadPosSingle(
		cnt4_t a0, cnt4_t c0, cnt4_t g0, cnt4_t t0,
		cnt4_t a1, cnt4_t c1, cnt4_t g1, cnt4_t t1,
		cnt4_t a2, cnt4_t c2, cnt4_t g2, cnt4_t t2,
		cnt4_t a3, cnt4_t c3, cnt4_t g3, cnt4_t t3,
		array<score_t, 4> const& pi) noexcept {

		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt4_t const r0 = a0 + g0, y0 = c0 + t0;
		cnt4_t const r1 = a1 + g1, y1 = c1 + t1;
		cnt4_t const r2 = a2 + g2, y2 = c2 + t2;
		cnt4_t const r3 = a3 + g3, y3 = c3 + t3;

		cnt4_t const rryy = quadXXYY(r0, r1, r2, r3, y0, y1, y2, y3);

		cnt4_t const aayy = quadXXYY(a0, a1, a2, a3, y0, y1, y2, y3);
		cnt4_t const ggyy = quadXXYY(g0, g1, g2, g3, y0, y1, y2, y3);
		cnt4_t const rrcc = quadXXYY(r0, r1, r2, r3, c0, c1, c2, c3);
		cnt4_t const rrtt = quadXXYY(r0, r1, r2, r3, t0, t1, t2, t3);

		cnt4_t const aacc = quadXXYY(a0, a1, a2, a3, c0, c1, c2, c3);
		cnt4_t const aatt = quadXXYY(a0, a1, a2, a3, t0, t1, t2, t3);
		cnt4_t const ggcc = quadXXYY(g0, g1, g2, g3, c0, c1, c2, c3);
		cnt4_t const ggtt = quadXXYY(g0, g1, g2, g3, t0, t1, t2, t3);

		return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
			+ (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
	}

	inline static score_t quadPos(array<cnt_t, 4> const& cnt0, array<cnt_t, 4> const& cnt1,
		array<cnt_t, 4> const& cnt2, array<cnt_t, 4> const& cnt3, array<score_t, 4> const& pi) noexcept {

		return quadPosSingle(
			(cnt4_t)cnt0[0], (cnt4_t)cnt0[1], (cnt4_t)cnt0[2], (cnt4_t)cnt0[3],
			(cnt4_t)cnt1[0], (cnt4_t)cnt1[1], (cnt4_t)cnt1[2], (cnt4_t)cnt1[3],
			(cnt4_t)cnt2[0], (cnt4_t)cnt2[1], (cnt4_t)cnt2[2], (cnt4_t)cnt2[3],
			(cnt4_t)cnt3[0], (cnt4_t)cnt3[1], (cnt4_t)cnt3[2], (cnt4_t)cnt3[3],
			pi);
	}

	inline static array<score_t, 3> quadPos(array<array<cnt_t, 4>, 4> const& cnt, array<score_t, 4> const& pi) noexcept {
		return { quadPos(cnt[0], cnt[1], cnt[2], cnt[3], pi),
				quadPos(cnt[0], cnt[2], cnt[1], cnt[3], pi),
				quadPos(cnt[0], cnt[3], cnt[1], cnt[2], pi) };
	}

	// +++ CASTER_TRI: quadPos with R-subtraction (only count quartets containing R) +++
	inline static array<score_t, 3> quadPosMinusR(array<array<cnt_t, 4>, 4> const& cnt, array<array<cnt_t, 4>, 4> const& rCnt, array<score_t, 4> const& pi) noexcept {
		array<score_t, 3> all = quadPos(cnt, pi);

		auto g = [&](int c, int n) -> cnt4_t { return (cnt4_t)cnt[c][n] - (cnt4_t)rCnt[c][n]; };

		score_t minusR0 = quadPosSingle(g(0,0),g(0,1),g(0,2),g(0,3), g(1,0),g(1,1),g(1,2),g(1,3), g(2,0),g(2,1),g(2,2),g(2,3), g(3,0),g(3,1),g(3,2),g(3,3), pi);
		score_t minusR1 = quadPosSingle(g(0,0),g(0,1),g(0,2),g(0,3), g(2,0),g(2,1),g(2,2),g(2,3), g(1,0),g(1,1),g(1,2),g(1,3), g(3,0),g(3,1),g(3,2),g(3,3), pi);
		score_t minusR2 = quadPosSingle(g(0,0),g(0,1),g(0,2),g(0,3), g(3,0),g(3,1),g(3,2),g(3,3), g(1,0),g(1,1),g(1,2),g(1,3), g(2,0),g(2,1),g(2,2),g(2,3), pi);

		return { all[0] - minusR0, all[1] - minusR1, all[2] - minusR2 };
	}
	// +++ END CASTER_TRI +++

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}
	
	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}
	
	// +++ CASTER_TRI: elementScore with R-subtraction (each position's rColorCnts is for the correct ref) +++
	score_t elementScore(size_t iElement) const noexcept{
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		score_t res = 0;
		for (index_t iPos : iota((index_t)0, nPos)){
			index_t gPos = iGenomePosBegin + iPos;
			res += scorePosMinusR(colorCnts[gPos], rColorCnts[gPos], element.eqFreqs);
		}
		return res;
	}

	// +++ CASTER_TRI: elementQuadripartitionScores with R-subtraction +++
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		array<score_t, 3> res = {0, 0, 0};
		for (index_t iPos : iota((index_t)0, nPos)) {
			index_t gPos = iGenomePosBegin + iPos;
			array<score_t, 3> part = quadPosMinusR(colorCnts[gPos], rColorCnts[gPos], element.eqFreqs);
			for (index_t i : iota((index_t)0, (index_t)3)) res[i] += part[i];
		}
		return res;
	}
	// +++ END CASTER_TRI +++

	Color(SharedConstData const& data) noexcept : sharedConstData(data), colorCnts(data.nGenomePos), rColorCnts(data.nGenomePos) {}

	template<typename DataClasses> friend DataClasses DriverHelper::read();
};

ChangeLog logDriverHelper("DriverHelper",
	"2026-02-04", "Chao Zhang", "Little code refactoring, no functional change", "patch",
	// +++ CASTER_TRI +++
	"2026-06-12", "Zuizhi Chen", "CASTER_TRI: add fasta2ref parsing, multi-fasta input, and multi-ref support", "minor");
	// +++ END CASTER_TRI +++

namespace DriverHelper {

using namespace std;

template<typename T, typename T2> array<T, 4>& operator+=(array<T, 4>& a, const array<T2, 4>& b) {
	for (int j = 0; j < 4; j++) {
		a[j] += b[j];
	}
	return a;
}

template<typename T> T sum(const array<T, 4>& cnt) {
	T result = 0;
	for (int j = 0; j < 4; j++) {
		result += cnt[j];
	}
	return result;
}

// +++ CASTER_TRI: rewritten read() with multi-ref fasta2ref support +++
// Parses a fasta2ref file (each line: <fasta_path> <reference_species>),
// reads all fasta files, and aggregates into a single SharedConstData.
// Each element carries its own reference species; per-ref counts tracked for R-subtraction.
template<typename DataClass> DataClass read() {
	using cnt_taxon_t = DataClass::ParentClass::cnt_taxon_t;
	using cnt_t = DataClass::ParentClass::cnt_t;

	common::LogInfo log(1);
	log.log() << "Parsing fasta2ref file and reading input..." << std::endl;
	DataClass sharedConstData;

	// +++ CASTER_TRI: parse fasta2ref file +++
	string fasta2refFile = ARG.get<string>("input");
	ifstream finF2R(fasta2refFile);
	if (!finF2R.is_open()) {
		throw std::logic_error("Cannot open fasta2ref file: " + fasta2refFile);
	}

	// Get the directory of fasta2ref file for resolving relative paths
	string fasta2refDir;
	{
		size_t pos = fasta2refFile.find_last_of("/\\");
		if (pos != string::npos) fasta2refDir = fasta2refFile.substr(0, pos + 1);
	}

	// +++ CASTER_TRI: collect (fasta_path, ref_taxon_id) pairs; dedup refs for priorityTaxa +++
	vector<string> fastaFiles;
	vector<size_t> fileRefTaxonIds;
	vector<size_t>& priorityTaxa = sharedConstData.priorityTaxa;
	unordered_set<size_t> seenRefs;

	string fastaPath, refSpeciesName;
	while (finF2R >> fastaPath >> refSpeciesName) {
		if (fastaPath[0] != '/' && fastaPath[0] != '\\' && (fastaPath.size() < 2 || fastaPath[1] != ':')) {
			fastaPath = fasta2refDir + fastaPath;
		}
		fastaFiles.push_back(fastaPath);
		size_t refId = common::taxonName2ID[refSpeciesName];
		fileRefTaxonIds.push_back(refId);
		if (!seenRefs.count(refId)) {
			seenRefs.insert(refId);
			priorityTaxa.push_back(refId);
		}
	}
	finF2R.close();

	if (fastaFiles.empty()) {
		throw std::logic_error("No fasta files found in fasta2ref file: " + fasta2refFile);
	}

	log.log() << "#Reference species (distinct): " << priorityTaxa.size() << std::endl;
	for (size_t rId : priorityTaxa) log.log() << "  " << common::taxonName2ID[rId] << std::endl;
	log.log() << "#Fasta files: " << fastaFiles.size() << std::endl;
	// +++ END CASTER_TRI +++

	// +++ CASTER_TRI: create temporary fasta list file for AlignmentParser +++
	string tempListFile = "/tmp/caster_tri_fasta_list_" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())) + ".txt";
	{
		ofstream fout(tempListFile);
		for (const string& f : fastaFiles) fout << f << "\n";
	}
	// +++ END CASTER_TRI +++

	// +++ CASTER_TRI: iterate alignments; track ref for each file +++
	aligment_utilities::AlignmentParser AP(tempListFile, 2), AP2(tempListFile, 3);
	size_t iFile = 0;
	// +++ END CASTER_TRI +++

	while (AP.nextAlignment()) {
		// +++ CASTER_TRI: get the reference species for this alignment +++
		size_t fileRefTaxonId = fileRefTaxonIds[iFile++];
		// +++ END CASTER_TRI +++

		size_t nSites = AP.getLength();
		size_t chunkMaxSize = ARG.get<size_t>("chunk");
		size_t nChunk = (nSites + chunkMaxSize - 1) / chunkMaxSize;
		vector<vector<size_t> > sites(nChunk);
		vector<array<double, 4> > eqfreq;
		size_t iElementBegin = sharedConstData.elements.size();
		unordered_map<size_t, size_t> taxon2row;
			
		{
			size_t nTotalSpeciesmen = 0;
			unordered_map<size_t, size_t> nSpeciesmen;
			vector<array<unsigned short, 4> > freq;
			freq.resize(AP.getLength());
			while (AP.nextSeq()) {
				size_t iTaxon = common::taxonName2ID[AP.getName()];
				nTotalSpeciesmen++;
				nSpeciesmen[iTaxon]++;

				if (!taxon2row.count(iTaxon)) taxon2row[iTaxon] = taxon2row.size();
				string seq = AP.getSeq();
				for (size_t i = 0; i < seq.size(); i++) {
					switch (seq[i]) {
						case 'A': freq[i][0]++; break;
						case 'C': freq[i][1]++; break;
						case 'G': freq[i][2]++; break;
						case 'T': freq[i][3]++; break;
					}
				}
			}

			size_t maxSpeciesman = 0;
			for (auto const& element : nSpeciesmen) {
				maxSpeciesman = std::max(maxSpeciesman, element.second);
			}

			if (std::same_as<cnt_taxon_t, bool> && maxSpeciesman >= 2) {
				log.log() << "Seems there is more than one haploid genome per taxon and thus bool type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_taxon_t, unsigned char> && maxSpeciesman >= 256) {
				log.log() << "Seems there are more than 255 haploid genomes per taxon (which is fishy) and thus unsigned char type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_taxon_t, unsigned short> && maxSpeciesman >= 65536) {
				common::LogInfo err(-100);
				err.log() << "Seems there are more than 65535 haploid genomes per taxon (which is astonishing)! Please ask the author for a specially made version..." << std::endl;
				exit(-1);
			}
			if (std::same_as<cnt_t, unsigned char> && nTotalSpeciesmen >= 256) {
				log.log() << "Seems there are more than 255 haploid genomes in total and thus unsigned char type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_t, unsigned short> && nTotalSpeciesmen >= 65536) {
				common::LogInfo err(-100);
				err.log() << "Seems there are more than 65535 haploid genomes in total (which is fishy)! Please ask the author for a specially made version..." << std::endl;
				exit(-1);
			}

			for (size_t i = 0; i < nChunk; i++) {
				size_t s = i * nSites / nChunk, t = (i + 1) * nSites / nChunk;
				array<size_t, 4> sumFreq = {};
				for (size_t j = s; j < t; j++) {
					sumFreq += freq[j];
				#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
					sites[i].push_back(j);
				#else
					if (freq[j][0] + freq[j][2] >= 2 && freq[j][1] + freq[j][3] >= 2) sites[i].push_back(j);
				#endif
				}
				double total = sum(sumFreq);
				if (total > 0) eqfreq.push_back({ sumFreq[0] / total, sumFreq[1] / total, sumFreq[2] / total, sumFreq[3] / total });
				else eqfreq.push_back({ 0.25, 0.25, 0.25, 0.25 });
			}
		}

		AP2.nextAlignment();
		for (size_t i = 0; i < nChunk; i++) {
			typename DataClass::Element element;
			element.iGenomePosBegin = sharedConstData.nGenomePos;
			element.nPos = sites[i].size();
			element.cnts.resize(taxon2row.size(), vector<array<typename DataClass::ParentClass::cnt_taxon_t, 4> >(element.nPos));
			element.taxon2row.resize(common::taxonName2ID.nTaxa(), -1);
			element.eqFreqs = eqfreq[i];
			// +++ CASTER_TRI: tag element with its reference species +++
			element.iReferenceTaxonId = fileRefTaxonId;
			// +++ END CASTER_TRI +++
			sharedConstData.elements.push_back(element);
			sharedConstData.nGenomePos += element.nPos;
		}

		while (AP2.nextSeq()) {
			size_t iTaxon = common::taxonName2ID[AP2.getName()];
			size_t iRow = taxon2row[iTaxon];
			string seq = AP2.getSeq();
			for (size_t iChunk : iota((size_t) 0, nChunk)) {
				typename DataClass::Element &element = sharedConstData.elements[iElementBegin + iChunk];
				element.taxon2row[iTaxon] = iRow;
				for (size_t iPos : iota((size_t) 0, sites[iChunk].size())) {
					switch (seq[sites[iChunk][iPos]]) {
						case 'A': element.cnts[iRow][iPos][0] += 1; break;
						case 'C': element.cnts[iRow][iPos][1] += 1; break;
						case 'G': element.cnts[iRow][iPos][2] += 1; break;
						case 'T': element.cnts[iRow][iPos][3] += 1; break;
					}
				}
			}
		}
	}

	// +++ CASTER_TRI: clean up temp file +++
	std::remove(tempListFile.c_str());
	// +++ END CASTER_TRI +++

	// +++ CASTER_TRI: cache priority taxa for taxonOrderPrioritizing static method +++
	s_priorityTaxa = sharedConstData.priorityTaxa;
	// +++ END CASTER_TRI +++

	return sharedConstData;
}
// +++ END CASTER_TRI +++

};

ChangeLog logDriver("Driver",
	"2026-02-01", "Chao Zhang", "Change prgramName to caster", "patch",
	"2026-02-08", "Chao Zhang", "Adding more type support", "patch",
	// +++ CASTER_TRI: update driver for caster-tri +++
	"2026-06-12", "Zuizhi Chen", "CASTER_TRI: rename to caster-tri, use -i for fasta2ref input", "minor");
	// +++ END CASTER_TRI +++

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<bool, unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char, unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<bool, unsigned short> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char, unsigned short> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned short, unsigned short> >::SharedConstData>;

	// +++ CASTER_TRI: updated program name +++
	static std::pair<string, string> programNames() noexcept {
		return { "caster-tri", "CASTER-TRI: Coalescence-aware Alignment-based Species Tree EstimatoR with TRI-reference" };
	}
	// +++ END CASTER_TRI +++

	static void addArguments() noexcept {
		ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		try { return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<1, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<2, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<3, DataClasses> >(); } catch (...) {}
		return DriverHelper::read<std::variant_alternative_t<4, DataClasses> >();
	}
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# CASTER-TRI: Coalescence-aware Alignment-based Species Tree EstimatoR with TRI-reference

[<img src="../misc/CASTER.png" width="500"/>](../misc/CASTER.png)

CASTER-TRI is a variant of CASTER that restricts quartet scoring to only those quartets that contain a triangulation reference species (R). This is designed for phylogenomic analyses where multiple fasta files are aligned to a common reference genome.

## How it works
1. Multiple fasta alignment files are provided via a `fasta2ref` file.
2. A single reference species (R) is specified per alignment batch.
3. Only quartets containing the reference species contribute to the tree score.
4. The reference species is always placed first during heuristic search to avoid information loss.

## Key difference from CASTER
While CASTER scores all possible quartets across the alignment, CASTER-TRI only scores quartets that include the designated reference species. This makes it suitable for analyses where the reference genome provides a consistent triangulation anchor.

This is a modification from CASTER (Coalescence-aware Alignment-based Species Tree EstimatoR). See the CASTER publication for details of the underlying statistical model.

## Publication

Chao Zhang, Rasmus Nielsen, Siavash Mirarab, CASTER: Direct species tree inference from whole-genome alignments. Science (2025) https://www.science.org/doi/10.1126/science.adk9688
)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT

CASTER-TRI uses a `fasta2ref` file to specify the mapping between fasta alignment files and their reference species.

## fasta2ref file format
Each line contains a fasta file path followed by the reference species name, separated by whitespace:
```
/path/to/alignment1.fasta    Felis_catus
/path/to/alignment2.fasta    Felis_catus
```
Currently, all alignment files must use the same reference species (single-reference mode).

## Fasta file format
* Each fasta file should be a multiple sequence alignment in FASTA format.
* All fasta files should contain the same set of species (may be in different order).
* Each fasta file represents the same genomic region aligned to a different reference genome or using different alignment parameters.
* Only A, C, G, T characters are considered. Gap characters ('-') are ignored.

Example fasta file:
```
>Felis_catus
AACCTTGG
>Panthera_tigris
AACCTTGG
>Panthera_pardus
AACCGTGG
>Acinonyx_jubatus
AACCGTCG
```

Multiple speciesmen per species are supported for multiploid data:
```
>Felis_catus
AACCTTGG
>Felis_catus
AACCTTGC
```

## Usage
```
bin/caster-tri --fasta2ref example/fasta2ref.txt -o output.tre
```
)YOHANETYO";
	}

	// +++ CASTER_TRI +++
	string programName() const noexcept override { return "caster-tri"; }
	// +++ END CASTER_TRI +++

	string exampleInput() const noexcept override { return "fasta2ref.txt"; }
};

};
#endif
