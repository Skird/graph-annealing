#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <algorithm>

#include <boost/algorithm/string.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

namespace po = boost::program_options;

typedef float float_t;

struct Point {
	float_t x, y;

	Point() {}
	Point(float_t nx, float_t ny): x(nx), y(ny) {}

	Point operator - (const Point &p) const {
		return Point(x - p.x, y - p.y);
	}

	Point operator + (const Point &p) const {
		return Point(x + p.x, y + p.y);
	}

	float_t operator * (const Point &p) const {
		return x * p.x + y * p.y;
	}

	Point operator * (const float_t k) const {
		return Point(x * k, y * k);
	}
};

struct GraphLayout {
	int n;
	std::vector< std::vector<int> > adj;
	std::vector<Point> layout;
	float_t W, H;

	GraphLayout(
			int n,
			float_t nodeDistW, float_t borderW, float_t edgeW,
			float_t crossW, float_t crossT,
			float_t canvasW, float_t canvasH
	) {
		this->n = n;
		_nw = nodeDistW;
		_bw = borderW;
		_ew = edgeW;
		_cw = crossW;
		_mind = crossT;
		W = canvasW;
		H = canvasH;
		adj.assign(n, std::vector<int>());
	}

	void addEdge(int u, int v) {
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	template<typename RNG> void initialize_positions(const std::vector<Point> &start, RNG &rng) {
		layout = start;
		if (layout.empty()) {
			layout.resize(n);
			std::uniform_real_distribution<float_t> wdis(0, W);
			std::uniform_real_distribution<float_t> hdis(0, H);
			for (int i = 0; i < n; ++i) {
				layout[i].x = wdis(rng);
				layout[i].y = hdis(rng);
			}
		}
		initEnergy();
	}

	void move(int v, Point d) {
		_tmpne = _ne;
		_tmpbe = _be;
		_tmpee = _ee;
		_tmpce = _ce;
		_tmpPoint = layout[v];
		addEnergy(v, -1);
		layout[v].x = std::clamp<float_t>(
				layout[v].x + d.x, W / 1000, static_cast<const float_t &>(W * 0.999)
		);
		layout[v].y = std::clamp<float_t>(
				layout[v].y + d.y, H / 1000, static_cast<const float_t &>(H * 0.999)
		);
		addEnergy(v, 1);
	}

	// Take care to call this only straight after move
	void undoMove(int v, Point d) {
		layout[v] = _tmpPoint;
		_ne = _tmpne;
		_be = _tmpbe;
		_ee = _tmpee;
		_ce = _tmpce;
	}

	float_t getEnergy() const {
		return _ne + _be + _ee + _ce;
	}

	float_t getNodeDistEnergy() const {
		return _ne;
	}

	float_t getBorderEnergy() const {
		return _be;
	}

	float_t getEdgesEnergy() const {
		return _ee;
	}

	float_t getCrossingsEnergy() const {
		return _ce;
	};

private:
	float_t _mind;
	float_t _ne, _be, _ee, _ce;
	float_t _tmpne, _tmpbe, _tmpee, _tmpce;
	Point _tmpPoint;
	float_t _nw, _bw, _ew, _cw;

	float_t _sqr(float_t x) {
		return x * x;
	}

	float_t _sqdist(int i, int j) {
		return _sqr(layout[i].x - layout[j].x) + _sqr(layout[i].y - layout[j].y);
	}

	float_t _segsqdist(int u, int v, int w) {
		float_t ans = fmin(_sqdist(w, u), _sqdist(w, v));
		Point uw = layout[w] - layout[u];
		Point uv = layout[v] - layout[u];
		float_t cp = uw * uv;
		if (cp < 0) return ans;
		if ((layout[u] - layout[v]) * (layout[w] - layout[v]) < 0) return ans;
		cp /= uv * uv;
		Point h = layout[u] + (uv * cp);
		return _sqr(layout[w].x - h.x) + _sqr(layout[w].y - h.y);
	}

	void initEnergy() {
		_ne = _be = _ee = _ce = 0;
		if (_nw != 0) {
			float_t dE = 0;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < i; ++j) {
					dE += 1 / _sqdist(i, j);
				}
			}
			_ne += dE * _nw;
		}
		if (_bw != 0) {
			float_t dE = 0;
			for (int i = 0; i < n; ++i) {
				dE +=
						1 / _sqr(layout[i].x) + 1 / _sqr(W - layout[i].x) +
						1 / _sqr(layout[i].y) + 1 / _sqr(H - layout[i].y);
			}
			_be += _bw * dE;
		}
		if (_ew != 0) {
			float_t dE = 0;
			for (int v = 0; v < n; ++v) {
				for (const auto &u : adj[v]) {
					if (v < u) {
						dE += _sqdist(v, u);
					}
				}
			}
			_ee += _ew * dE;
		}
		if (_cw != 0) {
			float_t dE = 0;
			for (int v = 0; v < n; ++v) {
				for (int u : adj[v]) {
					if (u < v) continue;
					for (int w = 0; w < n; ++w) {
						if (v == w || u == w) continue;
						dE += 1 / fmax(_mind, _segsqdist(u, v, w));
					}
				}
			}
			_ce += _cw * dE;
		}
	}

	void addEnergy(int v, double c) {
		if (_nw != 0) {  // exact equality is fine here
			float_t dE = 0;
			for (int i = 0; i < n; ++i) {
				if (i != v) dE += 1 / _sqdist(v, i);
			}
			_ne += c * _nw * dE;
		}
		if (_bw != 0) {
			float_t dE =
					1 / _sqr(layout[v].x) + 1 / _sqr(W - layout[v].x) +
					1 / _sqr(layout[v].y) + 1 / _sqr(H - layout[v].y);
			_be += c * _cw * dE;
		}
		if (_ew != 0) {
			float_t dE = 0;
			for (const auto &u: adj[v]) {
				dE += _sqdist(v, u);
			}
			_ee += c * _ew * dE;
		}
		if (_cw != 0) {
			float_t dE = 0;
			for (int u = 0; u < n; ++u) {
				for (int w : adj[u]) {
					if (u == v || w == v || w < u) continue;
					dE += 1 / fmax(_mind, _segsqdist(u, w, v));
				}
			}
			for (int u : adj[v]) {
				for (int w = 0; w < n; ++w) {
					if (v == w || u == w) continue;
					dE += 1 / fmax(_mind, _segsqdist(v, u, w));
//					std::cerr << "DEBUG " << _segsqdist(v, u, w) << "\n";
//					std::cerr << "DEBUG " << _segsqdist(u, v, w) << "\n";
//					std::cerr << "====\n";
				}
			}
			_ce += c * _cw * dE;
		}
	}
};

float_t transition_probability(float_t dE, float_t T) {
	// here dE must be positive
	return exp(-dE / T);
}

void print_layout(const GraphLayout &layout) {
	for (int i = 0; i < layout.n; ++i) {
		std::cout << layout.layout[i].x << " " << layout.layout[i].y << ";";
	}
	std::cout
		<< layout.getNodeDistEnergy() << ";" << layout.getBorderEnergy()
		<< ";" << layout.getEdgesEnergy() << ";" << layout.getCrossingsEnergy() << "\n";
}

template<typename RNG>
void simulate_annealing(
		GraphLayout &layout, int iters, float_t startT,
		float_t dT, bool fine_tune, RNG &rng, bool verbose
) {
	std::uniform_int_distribution<int> vdis(0, layout.n - 1);
	std::uniform_real_distribution<float_t> pdis(0, 1);

	auto curS = static_cast<float_t>(hypotl(layout.W, layout.H) / 2);
	auto dS = static_cast<float_t>(powl(0.001, (float_t) 1.0 / iters));
	float_t temp = startT;

	if (fine_tune) {
		curS = static_cast<float_t>(hypotl(layout.W, layout.H) / 300);
		dS = 1.0;
	}

	print_layout(layout);
	for (int i = 0; i < iters; ++i) {
		int v = vdis(rng);
		Point d{};
		d.x = (2 * std::generate_canonical<float_t, 1, RNG>(rng) - 1) * curS;
		d.y = (2 * std::generate_canonical<float_t, 1, RNG>(rng) - 1) * curS;
		float_t was_energy = layout.getEnergy();
		if (verbose) {
			std::cerr << "Iteration #" << i << ": E=" << was_energy << " T=" << temp << std::endl;
			std::cerr << "ne=" << layout.getNodeDistEnergy();
			std::cerr << " be=" << layout.getBorderEnergy();
			std::cerr << " ee=" << layout.getEdgesEnergy();
			std::cerr << " ce=" << layout.getCrossingsEnergy() << std::endl;
		}
		layout.move(v, d);
		float_t new_energy = layout.getEnergy();
		if (
				(fine_tune && new_energy > was_energy) ||
				(new_energy > was_energy &&
					pdis(rng) > transition_probability(new_energy - was_energy, temp)
				)
		) {
			layout.undoMove(v, d);
			if (verbose) {
				std::cerr
					<< "undo occured (dE=" << new_energy - was_energy
					<< ", T=" << temp << ") prob="
					<< 1 - transition_probability(new_energy - was_energy, temp) << ")\n";
			}
		}
		temp *= dT;
		curS *= dS;
		print_layout(layout);
	}
}

int main(int argc, char **argv) {
	freopen("graph.txt", "r", stdin);
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "show usage")
		("nodes", po::value<float_t>(), "weight of even node distribution in energy")
		("edges", po::value<float_t>(), "weight of length of edges in energy")
		("border", po::value<float_t>(), "weight of proximity to borderline in energy")
		("cross", po::value<float_t>(), "weight of crossings in energy")
		("cross-threshold", po::value<float_t>(), "minimum distance value for node edge interaction")
		("startT", po::value<float_t>(), "statring temperature")
		("cooling", po::value<float_t>(), "schedule of temperature reduction")
		("rounds", po::value<int>(), "number of rounds of annealing")
		("width", po::value<float_t>(), "width of canvas")
		("height", po::value<float_t>(), "height of canvas")
		("starting-layout", po::value<std::string>(), "path to starting layout file (if omitted, random one is used)")
		("seed", po::value<int>(), "random seed")
		("fine-tune", "flag to enable local search (used for fine tuning)")
		("verbose", "verbose mode");
	po::variables_map varmap;
	po::store(po::parse_command_line(argc, argv, desc), varmap);
	po::notify(varmap);

	if (varmap.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}

	bool failed = false;
	for (const auto &s: {"nodes", "edges", "border", "startT", "cooling", "rounds", "width", "height"}) {
		if (!varmap.count(s)) {
			std::cerr << "Option --" << s << " is required" << std::endl;
			failed = true;
		}
	}
	if (failed) {
		std::cout << desc << std::endl;
		return 1;
	}

	int n, m;
	std::cin >> n >> m;
	GraphLayout layout(
			n,
			varmap["nodes"].as<float_t>(),
			varmap["border"].as<float_t>(),
			varmap["edges"].as<float_t>(),
			(varmap.count("cross") ? varmap["cross"].as<float_t>() : 0),
			(varmap.count("cross-threshold") ? varmap["cross-threshold"].as<float_t>() : 1),
			varmap["width"].as<float_t>(),
			varmap["height"].as<float_t>()
	);
	for (int i = 0; i < m; ++i) {
		int from, to;
		std::cin >> from >> to;
		layout.addEdge(from, to);
	}

	std::vector<Point> start;
	if (varmap.count("starting-layout")) {
		std::ifstream start_file(varmap["starting-layout"].as<std::string>());
		std::string s, last;
		while (std::getline(start_file, s)) {
			last = s;
		}
		std::istringstream ss(last);
		start.resize(static_cast<unsigned long>(n));
		for (int i = 0; i < n; ++i) {
			char sep;
			ss >> start[i].x >> start[i].y >> sep;
		}
	}
	std::mt19937 rng(varmap.count("seed") ? varmap["seed"].as<int>() : 0);
	layout.initialize_positions(start, rng);

	simulate_annealing(
			layout, varmap["rounds"].as<int>(),
			varmap["startT"].as<float_t>(), varmap["cooling"].as<float_t>(),
			static_cast<bool>(varmap.count("fine-tune")), rng,
			static_cast<bool>(varmap.count("verbose"))
	);
}
