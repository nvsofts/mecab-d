module mecab;

import std.stdio;
import std.array;
import std.string;
import std.conv;
import std.algorithm;
import core.stdc.config;

extern (C)
{
  alias mecab_t = void;
  alias mecab_lattice_t = void;

  struct mecab_path_t
  {
    mecab_node_t* rnode;
    mecab_path_t* rnext;
    mecab_node_t* lnode;
    mecab_path_t* lnext;
    int cost;
    float prob;
  }

  struct mecab_node_t
  {
    mecab_node_t* prev;
    mecab_node_t* next;
    mecab_node_t* enext;
    mecab_node_t* bnext;
    mecab_path_t* rpath;
    mecab_path_t* lpath;

    const char* surface;
    const char* feature;

    uint id;
    ushort length;
    ushort rlength;
    ushort rcAttr;
    ushort lcAttr;
    ushort posid;
    ubyte char_type;
    ubyte stat;
    ubyte isbest;

    float alpha;
    float beta;
    float prob;
    short wcost;
    c_long cost;
  }

  mecab_t* mecab_new(int argc, char** argv);
  mecab_t* mecab_new2(const char* arg);
  void mecab_destroy(mecab_t* mecab);

  char* mecab_version();

  char* mecab_strerror(mecab_t* mecab);

  int mecab_get_partial(mecab_t* mecab);
  void mecab_set_partial(mecab_t* mecab, int partial);

  float mecab_get_theta(mecab_t* mecab);
  void mecab_set_theta(mecab_t* mecab, float theta);

  int mecab_get_lattice_level(mecab_t* mecab);
  void mecab_set_lattice_level(mecab_t* mecab, int level);

  int mecab_get_all_morphs(mecab_t* mecab);
  void mecab_set_all_morphs(mecab_t* mecab, int all_morphs);

  int mecab_parse_lattice(mecab_t* mecab, mecab_lattice_t* lattice);

  char* mecab_sparse_tostr(mecab_t* mecab, const char* str);
  char* mecab_sparse_tostr2(mecab_t* mecab, const char* str, size_t len);
  char* mecab_sparse_tostr3(mecab_t* mecab, const char* str, size_t len, char* ostr, size_t olen);

  mecab_node_t* mecab_sparse_tonode(mecab_t *mecab, const char* str);

  // lattice
  mecab_lattice_t* mecab_lattice_new();
  void mecab_lattice_destroy(mecab_lattice_t *lattice);

  void mecab_lattice_clear(mecab_lattice_t *lattice);
}

class Node
{
  private mecab_node_t* node;

  this(mecab_node_t* node)
  {
    this.node = node;
  }

  @property mecab_node_t* ptr()
  {
    return this.node;
  }

  @property string surface()
  {
    return to!string(this.node.surface)[0..this.node.length];
  }

  @property string feature()
  {
    return to!string(this.node.feature);
  }

  struct Range
  {
    Node outer;

    this(mecab_node_t* node)
    {
      if (node is null) {
        this.outer = null;
        return;
      }

      this.outer = new Node(node);
    }

    @property bool empty()
    {
      return (this.outer is null);
    }

    @property Node front()
    {
      return this.outer;
    }

    void popFront()
    {
      if (this.outer.ptr.next is null) {
        this.outer = null;
        return;
      }

      this.outer = new Node(this.outer.ptr.next);
    }
  }
}

class MeCab
{
  private mecab_t* mecab;

  invariant
  {
    assert(this.mecab != null, "this.mecab is null!");
  }

  static string getVersion()
  {
    return to!string(mecab_version());
  }

  this(int argc, string[] argv)
  {
    const(char)*[] _argv;

    _argv = array(map!toStringz(argv));

    this.mecab = mecab_new(argc, cast(char**)_argv.ptr);

    if (this.mecab == null) {
      throw new Exception("Failed to initialize MeCab");
    }
  }

  this(string arg)
  {
    this.mecab = mecab_new2(toStringz(arg));

    if (this.mecab == null) {
      throw new Exception("Failed to initialize MeCab");
    }
  }

  ~this()
  {
    mecab_destroy(this.mecab);
  }

  @property mecab_t* ptr()
  {
    return this.mecab;
  }

  string strError()
  {
    return to!string(mecab_strerror(this.mecab));
  }

  @property int partial()
  {
    return mecab_get_partial(this.mecab);
  }

  @property int partial(int partial)
  {
    mecab_set_partial(this.mecab, partial);

    return partial;
  }

  @property float theta()
  {
    return mecab_get_theta(this.mecab);
  }

  @property float theta(float theta)
  {
    mecab_set_theta(this.mecab, theta);

    return theta;
  }

  @property int latticeLevel()
  {
    return mecab_get_lattice_level(this.mecab);
  }

  @property int latticeLevel(int level)
  {
    mecab_set_lattice_level(this.mecab, level);

    return level;
  }

  @property int allMorphs()
  {
    return mecab_get_all_morphs(this.mecab);
  }

  @property int allMorphs(int all_morphs)
  {
    mecab_set_all_morphs(this.mecab, all_morphs);

    return all_morphs;
  }

  int parseLattice(Lattice lattice)
  {
    return mecab_parse_lattice(this.mecab, lattice.ptr);
  }

  string sparseToStr(string str)
  {
    return to!string(mecab_sparse_tostr(this.mecab, toStringz(str)));
  }

  Node.Range sparseToNode(string str)
  {
    return Node.Range(mecab_sparse_tonode(this.mecab, toStringz(str)));
  }
}

class Lattice
{
  private mecab_lattice_t* lattice;

  invariant
  {
    assert(this.lattice != null, "this.lattice is null!");
  }

  this()
  {
    this.lattice = mecab_lattice_new();
  }

  ~this()
  {
    mecab_lattice_destroy(this.lattice);
  }

  @property mecab_lattice_t* ptr()
  {
    return this.lattice;
  }

  void clear()
  {
    mecab_lattice_clear(this.lattice);
  }
}

// for test
void main(string[] args)
{
  MeCab mecab = new MeCab(cast(int)args.length, args);

  foreach (n; mecab.sparseToNode("太郎は次郎が持っている本を花子に渡した。")) {
    writeln(n.surface);
    writeln(n.feature);
  }
}
