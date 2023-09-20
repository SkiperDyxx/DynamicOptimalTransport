#include "emd.hpp"
#include "pybind11/detail/common.h"
#include "pybind11/pytypes.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
#include <string>
namespace py = pybind11;

template <typename _FlowType, typename _CostType>
emd <py::array_t <_CostType>, _FlowType, _CostType, uint32_t> make_emd_py(py::array_t <_CostType> A, py::array_t <_FlowType> B, py::object dist) {
	if (A.ndim() != 2)
		throw std::runtime_error("Dimension of `point` array must be 2. Currently " + std::to_string(A.ndim()));
	if (B.ndim() != 1)
		throw std::runtime_error("Dimension of `flow` array must be 1. Currently " + std::to_string(B.ndim()));
	if (A.shape()[0] != B.shape()[0])
		throw std::runtime_error("Number of points mismatches. `point` array indicate " + std::to_string(A.shape()[0]) + " while `flow` array indicate " + std::to_string(B.shape()[0]));
	auto a = A.template unchecked<2>();
	auto b = B.template unchecked<1>();
	size_t n = a.shape(0), d = a.shape(1);
	std::vector <py::array_t <_CostType> > p; p.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		const _CostType* data = static_cast <const _CostType*>(a.data(i, 0));
		py::array_t <_CostType> vec(d);
		std::copy(data, data + d, vec.mutable_data());
		p.emplace_back(vec);
	}
#ifdef PYBIND11_DETAILED_ERROR_MESSAGES
	for (size_t i = 0; i < 3; ++i) {
		const _CostType *data = static_cast <const _CostType*>(p[i].data(0));
		for (size_t j = 0; j < 28; ++j) {
			for (size_t k = 0; k < 28; ++k)
				std::cout.put(data[j * 28 + k] ? '.' : ' ');
			std::cout.put('\n');
		}
		std::cout.put('\n');
	}
	std::cout.put('\n');
	for (size_t i = 0; i < 3; ++i) {
		const _CostType *data = static_cast <const _CostType*>(a.data(i, 0));
		for (size_t j = 0; j < 28; ++j) {
			for (size_t k = 0; k < 28; ++k)
				std::cout.put(data[j * 28 + k] ? '.' : ' ');
			std::cout.put('\n');
		}
		std::cout.put('\n');
	}
#endif
	const _FlowType *fdata = static_cast <const _FlowType*> (b.data(0));
	std::vector <_FlowType> f(fdata, fdata + n);
	return emd<py::array_t <_CostType>, _FlowType, _CostType, uint32_t>(std::move(p), std::move(f), [dist](py::array_t <_CostType> x, py::array_t <_CostType> y) -> _CostType { return dist(x, y).template cast<_CostType>(); });
}

PYBIND11_MODULE(dot, m) {
	m.doc() = "Python Dynamic Optimal Transport Library";
	using namespace py::literals;
	m.def("make_emd", &make_emd_py<int64_t, int64_t>, "Make Data Structure Instance", "point"_a, "flow"_a, "distfunc"_a);
	m.def("make_emd", &make_emd_py<int64_t, int32_t>, "Make Data Structure Instance", "point"_a, "flow"_a, "distfunc"_a);
	m.def("make_emd", &make_emd_py<int32_t, int64_t>, "Make Data Structure Instance", "point"_a, "flow"_a, "distfunc"_a);
	m.def("make_emd", &make_emd_py<int32_t, int32_t>, "Make Data Structure Instance", "point"_a, "flow"_a, "distfunc"_a);
	typedef emd <py::array_t <int64_t>, int64_t, int64_t, uint32_t> emd_64_64;
	typedef emd <py::array_t <int64_t>, int32_t, int64_t, uint32_t> emd_64_32;
	typedef emd <py::array_t <int32_t>, int64_t, int32_t, uint32_t> emd_32_64;
	typedef emd <py::array_t <int32_t>, int32_t, int32_t, uint32_t> emd_32_32;
	py::class_ <emd_64_64>(m, "emd_64_64")
		.def("cost", &emd_64_64::cost, "Get Current OT Cost")
		.def("modify", &emd_64_64::modify, "Move node in space", "idx"_a, "new_position"_a)
		.def("insert_demand", &emd_64_64::insert_demand, "Insert a node to demand set", "position"_a)
		.def("insert_supply", &emd_64_64::insert_supply, "Insert a node to supply set", "position"_a)
		.def("erase", &emd_64_64::erase, "Erase a point from node set", "idx"_a)
		.def("reserve", &emd_64_64::reserve, "Reserve Memory", "size"_a)
		.def("send_flow", &emd_64_64::send_flow, "Modify Weight", "from"_a, "to"_a, "Delta"_a);
	py::class_ <emd_64_32>(m, "emd_64_32")
		.def("cost", &emd_64_32::cost, "Get Current OT Cost")
		.def("modify", &emd_64_32::modify, "Move node in space", "idx"_a, "new_position"_a)
		.def("insert_demand", &emd_64_32::insert_demand, "Insert a node to demand set", "position"_a)
		.def("insert_supply", &emd_64_32::insert_supply, "Insert a node to supply set", "position"_a)
		.def("erase", &emd_64_32::erase, "Erase a point from node set", "idx"_a)
		.def("reserve", &emd_64_32::reserve, "Reserve Memory", "size"_a)
		.def("send_flow", &emd_64_32::send_flow, "Modify Weight", "from"_a, "to"_a, "Delta"_a);
	py::class_ <emd_32_64>(m, "emd_32_64")
		.def("cost", &emd_32_64::cost, "Get Current OT Cost")
		.def("modify", &emd_32_64::modify, "Move node in space", "idx"_a, "new_position"_a)
		.def("insert_demand", &emd_32_64::insert_demand, "Insert a node to demand set", "position"_a)
		.def("insert_supply", &emd_32_64::insert_supply, "Insert a node to supply set", "position"_a)
		.def("erase", &emd_32_64::erase, "Erase a point from node set", "idx"_a)
		.def("reserve", &emd_32_64::reserve, "Reserve Memory", "size"_a)
		.def("send_flow", &emd_32_64::send_flow, "Modify Weight", "from"_a, "to"_a, "Delta"_a);
	py::class_ <emd_32_32>(m, "emd_32_32")
		.def("cost", &emd_32_32::cost, "Get Current OT Cost")
		.def("modify", &emd_32_32::modify, "Move node in space", "idx"_a, "new_position"_a)
		.def("insert_demand", &emd_32_32::insert_demand, "Insert a node to demand set", "position"_a)
		.def("insert_supply", &emd_32_32::insert_supply, "Insert a node to supply set", "position"_a)
		.def("erase", &emd_32_32::erase, "Erase a point from node set", "idx"_a)
		.def("reserve", &emd_32_32::reserve, "Reserve Memory", "size"_a)
		.def("send_flow", &emd_32_32::send_flow, "Modify Weight", "from"_a, "to"_a, "Delta"_a);
}
