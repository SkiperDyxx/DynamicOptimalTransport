#pragma once
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include <tuple>
#include <list>

namespace lct {
	template <typename T = int64_t>
	struct node {
		node *fa, *son[2], *minp[2];
		T delta, minv[2];
		bool is_edge;
		bool edge_dir;
		// minv[1] 维护正向边
		// minv[0] 维护反向边
		// edge_dir = 1 表示正向边
		bool reverse_tag;
		// In push ups, the node need its childrens' minmax attributes
		// So minmax value should be changed once the reverse tag changes
		// Other attributes may change during push downs
	};
	template <typename T>
	static inline bool is_root(node<T> *p);
	template <typename T>
	static void access(node<T> *p);
	template <typename T>
	static void make_root(node<T> *p);
	template <typename T>
	static inline void push_up(node<T> *p);
	template <typename T>
	static inline void push_down(node<T> *p);
	template <typename T>
	static inline void add_val(node<T> *p, T val);
	template <typename T>
	static inline void reverse(node<T> *p);
	template <typename T>
	static void splay(node<T> *p);
	template <typename T>
	static inline void rotate(node<T> *p);
	template <typename T>
	node<T> *add_node();
	template <typename T>
	node<T> *link(node<T> *from, node<T> *to, T w = 0);
	template <typename T>
	void cut(node<T> *u, node<T> *v);
	template <typename T>
	std::pair <node<T>*, T> sendflow(node<T> *source, node<T> *target, T flow = std::numeric_limits<T>::max());

	template <typename T>
	node<T>* add_node() {
		node<T>* ptr = new node<T>;
		ptr->fa = ptr->son[0] = ptr->son[1] = nullptr;
		ptr->minp[0] = ptr->minp[1] = ptr;
		ptr->minv[0] = std::numeric_limits<T>::max();
		ptr->minv[1] = std::numeric_limits<T>::max();
		ptr->delta = 0;
		ptr->reverse_tag = false;
		ptr->edge_dir = false;
		ptr->is_edge = false;
		return ptr;
	}

	template <typename T>
	bool is_root(node<T>* p) {
		return p->fa == nullptr || (p->fa->son[0] != p && p->fa->son[1] != p);
	}

	template <typename T>
	void splay(node<T> *p) {
		static std::vector <node<T>*> stk;
		assert(stk.empty());
		for (node<T> *q = p; stk.push_back(q), !is_root(q); q = q->fa);
		while (!stk.empty()) {
			push_down(stk.back());
			stk.pop_back();
		}
		while (!is_root(p)) {
			node<T> *f = p->fa, *g = f->fa;
			if (!is_root(f))
				rotate((g->son[1] == f) != (f->son[1] == p) ? p : f);
			rotate(p);
		}
	}

	template <typename T>
	void add_val(node<T> *p, T val) {
		p->delta += val;
		if (p->minv[0] != std::numeric_limits<T>::max())
			p->minv[0] -= val;
		if (p->minv[1] != std::numeric_limits<T>::max())
			p->minv[1] += val;
	}

	template <typename T>
	void access(node<T> *p) {
		for (node<T> *t = nullptr; p != nullptr; t = p, p = p->fa) {
			splay(p);
			if (t != nullptr)
				add_val(t, T(-p->delta));
			if (p->son[1] != nullptr)
				add_val(p->son[1], p->delta);
			p->son[1] = t;
			push_up(p);
		}
	}

	template <typename T>
	void make_root(node<T> *p) {
		access(p);
		splay(p);
		assert(p->minv[0] >= 0);
		assert(p->minv[1] >= 0);
		reverse(p);
	}

	template <typename T>
	void reverse(node<T>* p) {
		p->delta = -p->delta;
		p->reverse_tag ^= 1;
		std::swap(p->minv[0], p->minv[1]);
		std::swap(p->minp[0], p->minp[1]);
	}

	template <typename T>
	void push_down(node<T> *p) {
		if (p->reverse_tag) {
			p->reverse_tag = 0;
			p->edge_dir ^= 1;
			std::swap(p->son[0], p->son[1]);
			for (int k = 0; k < 2; ++k)
				if (p->son[k] != nullptr) {
					reverse(p->son[k]);
					assert(p->son[k]->minv[0] == std::numeric_limits<T>::max() || p->minv[0] <= p->son[k]->minv[0] - p->delta);
					assert(p->son[k]->minv[1] == std::numeric_limits<T>::max() || p->minv[1] <= p->son[k]->minv[1] + p->delta);
				}
		}
	}

	template <typename T>
	void push_up(node<T> *p) {
		assert(!p->reverse_tag);
		p->minv[0] = std::numeric_limits<T>::max();
		p->minv[1] = std::numeric_limits<T>::max();
		if (p->is_edge) {
			if (p->edge_dir)
				p->minv[1] = p->delta;
			else
				p->minv[0] = -p->delta;
		}
		p->minp[0] = p->minp[1] = p;
		for (int k = 0; k < 2; ++k)
			if (p->son[k] != nullptr) {
				if (p->son[k]->minv[0] != std::numeric_limits<T>::max() && p->son[k]->minv[0] - p->delta < p->minv[0]) {
					p->minv[0] = p->son[k]->minv[0] - p->delta;
					p->minp[0] = p->son[k]->minp[0];
				}
				if (p->son[k]->minv[1] != std::numeric_limits<T>::max() && p->son[k]->minv[1] + p->delta < p->minv[1]) {
					p->minv[1] = p->son[k]->minv[1] + p->delta;
					p->minp[1] = p->son[k]->minp[1];
				}
			}
	}

	template <typename T>
	void rotate(node<T>* p) {
		node<T> *f = p->fa, *g = f->fa;
		bool k = f->son[1] == p;
		assert(f->son[k] == p);
		assert(p->reverse_tag == false);
		assert(f->reverse_tag == false);
		p->fa = g;
		if (g != nullptr) {
			if (g->son[0] == f)
				g->son[0] = p;
			else if (g->son[1] == f)
				g->son[1] = p;
		}
		if (p->son[!k] != nullptr) {
			p->son[!k]->fa = f;
			add_val(p->son[!k], p->delta);
		}
		f->son[k] = p->son[!k]; p->son[!k] = f; f->fa = p;
		p->delta += f->delta;
		f->delta -= p->delta;
		push_up(f); push_up(p);
	}

	template <typename T>
	node<T> *link(node<T> *u, node<T> *v, T w) {
		assert(w >= 0);
		make_root(v); access(u);
#ifndef NDEBUG
		for (node<T> *p = u; p != nullptr; p = p->fa)
			assert(p != v);
#endif
		node<T> *ptr = new node<T>;
		ptr->fa = u; ptr->son[0] = ptr->son[1] = nullptr;
		ptr->minv[1] = ptr->delta = w; ptr->minp[0] = ptr->minp[1] = ptr;
		ptr->minv[0] = std::numeric_limits<T>::max();
		ptr->reverse_tag = false;
		ptr->edge_dir = true;
		ptr->is_edge = true;
		v->fa = ptr;
		return ptr;
	}

	template <typename T>
	void cut(node<T> *u, node<T> *v) {
		make_root(u);
		access(v);
		splay(v);
		assert(v->son[1] == nullptr);
		assert(u->son[0] == nullptr);
		if (u->son[1] != nullptr) {
			assert(v->son[0] == u);
			assert(u->son[1]->son[0] == nullptr && u->son[1]->son[1] == nullptr);
			delete u->son[1];
			u->son[1] = u->fa = v->son[0] = nullptr;
			push_up(u); push_up(v);
		}
		else {
			assert(v->son[0]->son[0] == u);
			assert(v->son[0]->son[1] == nullptr);
			delete v->son[0];
			v->son[0] = u->fa = nullptr;
			push_up(v);
		}
	}

	template <typename T>
	std::pair <node<T>*, T> sendflow(node<T> *u, node<T> *v, T flow) {
		make_root(u);
		access(v);
		splay(v);
#ifndef NDEBUG
		{
			node<T> *p = u;
			while (!is_root(p))
				p = p->fa;
			assert(p == v);
		}
#endif
		assert(v->minv[0] >= 0);
		assert(v->minv[1] >= 0);
		flow = std::min(flow, v->minv[0]);
		add_val(v, flow);
		return std::make_pair(v->minp[0], flow);
	}

	template <typename T>
	T get_flow(node<T> *u, node<T> *v) {
		make_root(v);
		access(u);
		splay(u);
#ifndef NDEBUG
		{
			node<T> *p = v;
			while (!is_root(p))
				p = p->fa;
			assert(p == u);
		}
#endif
		assert(u->minv[0] >= 0);
		assert(u->minv[1] >= 0);
		return u->minv[0];
	}
}
