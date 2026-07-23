#ifndef WEIGHTED_TRIAL_HPP
#define WEIGHTED_TRIAL_HPP

#include "stepwise_colorable.hpp"
#include "alignment_utilities.hpp"

namespace weighted_trial{

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;
using std::unordered_map;

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

template<typename cnt_taxon_type = unsigned char> struct StepwiseColorDefaultAttributes {
	using score_t = double;
	using cnt_taxon_t = cnt_taxon_type;
	using cnt_t = double;
	using cnt4_t = double;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = 1e-3;
};

ChangeLog logColor("Color",
	"2026-07-06", "Zuizhi Chen", "Forked from CASTER_TRI v2.6.1", "minor",
	"2026-07-06", "Zuizhi Chen", "Sequence-similarity weighting with colorPairWeight for correct pair scoring", "minor");

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
			size_t iReferenceTaxonId = (size_t)-1; //reference species this element was aligned to
			vector<cnt_t> speciesWeights;        // WTRIAL: per-species-row sequence-similarity weight

			bool hasTaxon(size_t iTaxon) const noexcept{
				return iTaxon < taxon2row.size() && taxon2row[iTaxon] != -1;
			}
		};

		vector<Element> elements;
		index_t nGenomePos = 0;
		//all reference species taxon IDs used for priority ordering
		vector<size_t> priorityTaxa;


		size_t nElements() const noexcept { return elements.size(); }
	};

	//static method for TAXON_ORDER_PRIORITIZING concept
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

	//static cache of priority taxa, set during read()
	static inline std::vector<size_t> s_priorityTaxa;

private:
	SharedConstData const& sharedConstData;
	vector<array<array<cnt_t, 4>, 4>> colorWeight;     // WTRIAL: weighted-sum per pos/color/nuc
	vector<array<array<cnt_t, 4>, 4>> colorPairWeight; // WTRIAL: pair-weight sum per pos/color/nuc
	vector<array<cnt_t, 4>> refCnt;                    // refCnt[iGenomePos][iNucleotide] — reference species counts
	vector<size_t> refColor;                              // refColor[iGenomePos] — current color of ref, (size_t)-1 if not set

	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];
		if (!element.hasTaxon(iTaxon)) return;
		index_t iRow = element.taxon2row[iTaxon];
		index_t iPosBegin = element.iGenomePosBegin;

		if (iTaxon == element.iReferenceTaxonId) {
			for (index_t iPos : iota((index_t)0, element.nPos))
				refColor[iPosBegin + iPos] = (isSet) ? iColor : (size_t)-1;
			return;
		}

		cnt_t w = element.speciesWeights[iRow]; // WTRIAL: per-species similarity weight
		for (index_t iPos : iota((index_t)0, element.nPos)){
			index_t gPos = iPosBegin + iPos;
			for (index_t iNuc : iota((index_t)0, (index_t)4)) {
				cnt_taxon_t cnt = element.cnts[iRow][iPos][iNuc];
				if (cnt == 0) continue;
				cnt_t dw = (cnt_t)cnt * w;
				if constexpr (isSet) {
					// WTRIAL: add — order: pairWeight first (uses old colorWeight), then colorWeight
					colorPairWeight[gPos][iColor][iNuc] += colorWeight[gPos][iColor][iNuc] * dw;
					colorWeight[gPos][iColor][iNuc] += dw;
				} else {
					// WTRIAL: remove — order: colorWeight first, then pairWeight (uses new colorWeight)
					colorWeight[gPos][iColor][iNuc] -= dw;
					colorPairWeight[gPos][iColor][iNuc] -= colorWeight[gPos][iColor][iNuc] * dw;
				}
			}
		}
	}
	
	inline static cnt_t XXYY(cnt_t xR, cnt_t x0, cnt_t x1, cnt_t x2,
	                             cnt_t yR, cnt_t y0, cnt_t y1, cnt_t y2,
	                             cnt_t x11, cnt_t x22, cnt_t y11, cnt_t y22) noexcept{
		return xR * x0 * y1 * y2 * 2 + yR * y0 * x1 * x2 * 2 +
		       xR * (xR - 1) * y1 * y2 + yR * (yR - 1) * x1 * x2 +
		       xR * x1 * y22 * 2 + yR * y1 * x22 * 2 +
		       xR * x2 * y11 * 2 + yR * y2 * x11 * 2;
	}

	//WTRIAL: scorePos with pair-weighted formula
	inline static score_t scorePos(
		array<array<cnt_t, 4>, 4> const& cw,   // colorWeight (weighted sums)
		array<array<cnt_t, 4>, 4> const& cp,   // colorPairWeight (pair sums)
		array<cnt_t, 4> const& rCnt,           // ref counts
		array<score_t, 4> const& pi) noexcept
	{
		cnt_t const aR = rCnt[0], cR = rCnt[1], gR = rCnt[2], tR = rCnt[3];
		cnt_t const a0 = cw[0][0], c0 = cw[0][1], g0 = cw[0][2], t0 = cw[0][3];
		cnt_t const a1 = cw[1][0], c1 = cw[1][1], g1 = cw[1][2], t1 = cw[1][3];
		cnt_t const a2 = cw[2][0], c2 = cw[2][1], g2 = cw[2][2], t2 = cw[2][3];

		// WTRIAL: Pair-weights for purines (nuc 0,2) and pyrimidines (1,3) in each colour
		auto pw = [&](int col, int n1, int n2) -> cnt_t {
			if (n1 == n2) return cp[col][n1];
			return cp[col][n1] + cp[col][n2] + cw[col][n1] * cw[col][n2];
		};

		cnt_t const r11 = pw(1,0,2), r22 = pw(2,0,2);  // purine pairs (A+G cross-nuc)
		cnt_t const y11 = pw(1,1,3), y22 = pw(2,1,3);  // pyrimidine pairs (C+T cross-nuc)
		cnt_t const a11 = cp[1][0], a22 = cp[2][0];    // A pairs
		cnt_t const g11 = cp[1][2], g22 = cp[2][2];    // G pairs
		cnt_t const c11 = cp[1][1], c22 = cp[2][1];    // C pairs
		cnt_t const t11 = cp[1][3], t22 = cp[2][3];    // T pairs

		score_t const A = pi[0], Ci = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = Ci + T, R2 = A * A + G * G, Y2 = Ci * Ci + T * T;
		cnt_t const r0 = a0 + g0, y0 = c0 + t0;
		cnt_t const r1 = a1 + g1, y1 = c1 + t1;
		cnt_t const r2 = a2 + g2, y2 = c2 + t2;
		cnt_t const rR = aR + gR, yR = cR + tR;

		cnt_t const rryy = XXYY(rR, r0, r1, r2, yR, y0, y1, y2, r11, r22, y11, y22);
		cnt_t const aayy = XXYY(aR, a0, a1, a2, yR, y0, y1, y2, a11, a22, y11, y22);
		cnt_t const ggyy = XXYY(gR, g0, g1, g2, yR, y0, y1, y2, g11, g22, y11, y22);
		cnt_t const rrcc = XXYY(rR, r0, r1, r2, cR, c0, c1, c2, r11, r22, c11, c22);
		cnt_t const rrtt = XXYY(rR, r0, r1, r2, tR, t0, t1, t2, r11, r22, t11, t22);
		cnt_t const aacc = XXYY(aR, a0, a1, a2, cR, c0, c1, c2, a11, a22, c11, c22);
		cnt_t const aatt = XXYY(aR, a0, a1, a2, tR, t0, t1, t2, a11, a22, t11, t22);
		cnt_t const ggcc = XXYY(gR, g0, g1, g2, cR, c0, c1, c2, g11, g22, c11, c22);
		cnt_t const ggtt = XXYY(gR, g0, g1, g2, tR, t0, t1, t2, g11, g22, t11, t22);

		return (score_t)(rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
		     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y));
	}
	
	inline static cnt_t quadXXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t x3, cnt4_t y0, cnt4_t y1, cnt4_t y2, cnt4_t y3) noexcept {
		return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
	}

	// WTRIAL: cnt_t overloads for elementQuadripartitionScores
	static score_t quadPosSingle(array<cnt_t, 4> const& cnt0, array<cnt_t, 4> const& cnt1,
		array<cnt_t, 4> const& cnt2, array<cnt_t, 4> const& cnt3, array<score_t, 4> const& pi) noexcept {

		score_t const A = pi[0], Ci = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = Ci + T, R2 = A * A + G * G, Y2 = Ci * Ci + T * T;
		cnt_t const a0 = cnt0[0], c0 = cnt0[1], g0 = cnt0[2], t0 = cnt0[3], r0 = a0 + g0, y0 = c0 + t0;
		cnt_t const a1 = cnt1[0], c1 = cnt1[1], g1 = cnt1[2], t1 = cnt1[3], r1 = a1 + g1, y1 = c1 + t1;
		cnt_t const a2 = cnt2[0], c2 = cnt2[1], g2 = cnt2[2], t2 = cnt2[3], r2 = a2 + g2, y2 = c2 + t2;
		cnt_t const a3 = cnt3[0], c3 = cnt3[1], g3 = cnt3[2], t3 = cnt3[3], r3 = a3 + g3, y3 = c3 + t3;

		cnt_t const rryy = quadXXYY(r0, r1, r2, r3, y0, y1, y2, y3);
		cnt_t const aayy = quadXXYY(a0, a1, a2, a3, y0, y1, y2, y3);
		cnt_t const ggyy = quadXXYY(g0, g1, g2, g3, y0, y1, y2, y3);
		cnt_t const rrcc = quadXXYY(r0, r1, r2, r3, c0, c1, c2, c3);
		cnt_t const rrtt = quadXXYY(r0, r1, r2, r3, t0, t1, t2, t3);
		cnt_t const aacc = quadXXYY(a0, a1, a2, a3, c0, c1, c2, c3);
		cnt_t const aatt = quadXXYY(a0, a1, a2, a3, t0, t1, t2, t3);
		cnt_t const ggcc = quadXXYY(g0, g1, g2, g3, c0, c1, c2, c3);
		cnt_t const ggtt = quadXXYY(g0, g1, g2, g3, t0, t1, t2, t3);

		return (score_t)(rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
			+ (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y));
	}

	// WTRIAL: quadPos using cnt_t (double) arrays
	inline static array<score_t, 3> quadPos(array<array<cnt_t, 4>, 4> const& cnt, array<score_t, 4> const& pi) noexcept {
		return {quadPosSingle(cnt[0], cnt[1], cnt[2], cnt[3], pi),
		        quadPosSingle(cnt[0], cnt[2], cnt[1], cnt[3], pi),
		        quadPosSingle(cnt[0], cnt[3], cnt[1], cnt[2], pi)};
	}

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}
	
	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}
	
	//WTRIAL: Tripartition score with colorWeight/colorPairWeight (swap both arrays)
	score_t elementScore(size_t iElement) const noexcept{
		if (refColor[sharedConstData.elements[iElement].iGenomePosBegin] == (size_t)-1) return 0;
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		score_t res = 0;
		for (index_t iPos : iota((index_t)0, nPos)){
			index_t gPos = iGenomePosBegin + iPos;
			size_t C = refColor[gPos];
			if (C == (size_t)-1) continue;
			array<array<cnt_t, 4>, 4> cw = colorWeight[gPos];
			array<array<cnt_t, 4>, 4> cp = colorPairWeight[gPos];
			if (C != 0) { std::swap(cw[0], cw[C]); std::swap(cp[0], cp[C]); }
			res += scorePos(cw, cp, refCnt[gPos], element.eqFreqs);
		}
		return res;
	}

	//Edge NNI score — quadPos uses weighted sums directly (no pair-weight needed)
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		if (refColor[sharedConstData.elements[iElement].iGenomePosBegin] == (size_t)-1) return {0, 0, 0};
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		array<score_t, 3> res = {0, 0, 0};
		for (index_t iPos : iota((index_t)0, nPos)) {
			index_t gPos = iGenomePosBegin + iPos;
			size_t C = refColor[gPos];
			if (C == (size_t)-1) continue;
			array<array<cnt_t, 4>, 4> c = colorWeight[gPos];
			c[C] = refCnt[gPos];
			array<score_t, 3> part = quadPos(c, element.eqFreqs);
			for (index_t i : iota((index_t)0, (index_t)3)) res[i] += part[i];
		}
		return res;
	}

	//CASTER_TRI: constructor precomputes refCnt from elements;
	Color(SharedConstData const& data) noexcept : sharedConstData(data), colorWeight(data.nGenomePos), colorPairWeight(data.nGenomePos), refCnt(data.nGenomePos), refColor(data.nGenomePos, (size_t)-1) {
		for (auto const& element : sharedConstData.elements) {
			if (!element.hasTaxon(element.iReferenceTaxonId)) continue;
			index_t refRow = element.taxon2row[element.iReferenceTaxonId];
			for (index_t iPos : iota((index_t)0, element.nPos)) {
				index_t gPos = element.iGenomePosBegin + iPos;
				for (index_t nuc : iota((index_t)0, (index_t)4))
					refCnt[gPos][nuc] = element.cnts[refRow][iPos][nuc];
			}
		}
	}

	template<typename DataClasses> friend DataClasses DriverHelper::read();
};

ChangeLog logDriverHelper("DriverHelper",
	"2026-07-06", "Zuizhi Chen", "Forked from CASTER_TRI, compute per-species weights in read()", "minor");

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

//CASTER_TRI: fasta2ref parsing, multi-fasta input, multi-ref support
template<typename DataClass> DataClass read() {
	using cnt_taxon_t = DataClass::ParentClass::cnt_taxon_t;
	using cnt_t = DataClass::ParentClass::cnt_t;

	common::LogInfo log(1);
	log.log() << "Parsing fasta2ref file and reading input..." << std::endl;
	DataClass sharedConstData;

	string fasta2refFile = ARG.get<string>("input");
	ifstream finF2R(fasta2refFile);
	if (!finF2R.is_open()) {
		throw std::logic_error("Cannot open fasta2ref file: " + fasta2refFile);
	}

	string fasta2refDir;
	{
		size_t pos = fasta2refFile.find_last_of("/\\");
		if (pos != string::npos) fasta2refDir = fasta2refFile.substr(0, pos + 1);
	}

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

	//CASTER_TRI: create temporary fasta list file for AlignmentParser
	string tempListFile = "/tmp/caster_tri_fasta_list_" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())) + ".txt";
	{
		ofstream fout(tempListFile);
		for (const string& f : fastaFiles) fout << f << "\n";
	}

	aligment_utilities::AlignmentParser AP(tempListFile, 2), AP2(tempListFile, 3);
	size_t iFile = 0;

	while (AP.nextAlignment()) {
		size_t fileRefTaxonId = fileRefTaxonIds[iFile++];
		size_t nSites = AP.getLength();
		size_t chunkMaxSize = ARG.get<size_t>("chunk");
		size_t nChunk = (nSites + chunkMaxSize - 1) / chunkMaxSize;
		vector<vector<size_t> > sites(nChunk);
		vector<array<double, 4> > eqfreq;
		size_t iElementBegin = sharedConstData.elements.size();
		unordered_map<size_t, size_t> taxon2row;
			
		string refSeq; // WTRIAL: store ref sequence for weight computation
			
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
				if (iTaxon == fileRefTaxonId && refSeq.empty()) refSeq = seq; // WTRIAL: capture ref
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
			for (auto const& elem : nSpeciesmen) {
				maxSpeciesman = std::max(maxSpeciesman, elem.second);
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
			element.iReferenceTaxonId = fileRefTaxonId; //tag element with its reference species
			sharedConstData.elements.push_back(element);
			sharedConstData.nGenomePos += element.nPos;
		}

		// WTRIAL: compute per-species weights from Hamming distance to ref
		vector<cnt_t> speciesWeights(taxon2row.size(), 1.0);

		while (AP2.nextSeq()) {
			size_t iTaxon = common::taxonName2ID[AP2.getName()];
			size_t iRow = taxon2row[iTaxon];
			string seq = AP2.getSeq();

			// Weight: Hamming distance to ref (gaps excluded)
			cnt_t w = 1.0;
			if (iTaxon != fileRefTaxonId && !refSeq.empty()) {
				size_t hamming = 0, nonGap = 0, n = std::min(refSeq.size(), seq.size());
				for (size_t i = 0; i < n; i++) {
					if (refSeq[i] == '-' || seq[i] == '-') continue;
					nonGap++; if (refSeq[i] != seq[i]) hamming++;
				}
				double sim = (nonGap > 0) ? 1.0 - (double)hamming / (double)nonGap : 0.0;
				w = (sim < 0.25) ? 0.0 : (sim - 0.25) / 0.75;
			}
			speciesWeights[iRow] = w;
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
		// WTRIAL: assign weights to all elements of this alignment
		for (size_t iChunk = 0; iChunk < nChunk; iChunk++)
			sharedConstData.elements[iElementBegin + iChunk].speciesWeights = speciesWeights;
	}

	std::remove(tempListFile.c_str());

	DataClass::ParentClass::s_priorityTaxa = sharedConstData.priorityTaxa;

	return sharedConstData;
}

};

ChangeLog logDriver("Driver",
	"2026-07-06", "Zuizhi Chen", "Forked from CASTER_TRI v2.6.1", "minor");

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<bool> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned short> >::SharedConstData>;

	static std::pair<string, string> programNames() noexcept {
		return { "wtrial", "Weighted Trial: Weighted CASTER-TRI with sequence-similarity weighting" };
	}

	static void addArguments() noexcept {
		ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		try { return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<1, DataClasses> >(); } catch (...) {}
		return DriverHelper::read<std::variant_alternative_t<2, DataClasses> >();
	}
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# CASTER-TRI: Coalescence-aware Alignment-based Species Tree EstimatoR with TRI-reference

CASTER-TRI is a variant of CASTER that restricts quartet scoring to only those quartets that contain a reference species (R). This is designed for phylogenomic analyses where multiple fasta files are aligned to a common reference genome.

## How it works
1. Multiple fasta alignment files are provided via a `fasta2ref` file.
2. A reference species (R) is specified per alignment.
3. Only quartets containing the reference species contribute to the tree score.
4. The reference species is always placed first during heuristic search to avoid information loss.

## Key difference from CASTER
While CASTER scores all possible quartets across the alignment, CASTER-TRI only scores quartets that include the designated reference species. This makes it suitable for analyses where the reference genome provides a consistent triangulation anchor.

This is a modification from CASTER (Coalescence-aware Alignment-based Species Tree EstimatoR). See the CASTER publication for details of the underlying statistical model.

)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT

CASTER-TRI uses a `fasta2ref` file to specify the mapping between fasta alignment files and their reference species.

## fasta2ref file format
Each line contains a fasta file path followed by the reference species name, separated by whitespace:
```
/path/to/alignment1.fasta    Felis_catus
/path/to/alignment2.fasta    Otocolobus_manul
```

## Fasta file format
* Each fasta file should be a multiple sequence alignment in FASTA format.
* All fasta files should contain the same set of species (may be in different order).
* Each fasta file represents the same genomic region aligned to a different reference genome or using different alignment parameters.

## Usage
```
bin/caster-tri -i example/fasta2ref.txt -o output.tre
```
)YOHANETYO";
	}

	string programName() const noexcept override { return "wtrial"; }

	string exampleInput() const noexcept override { return "fasta2ref.txt"; }
};

};
#endif
