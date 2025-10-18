#ifndef LAHUTA_JSON_HPP
#define LAHUTA_JSON_HPP

#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gemmi/third_party/sajson.h"

// clang-format off
namespace lahuta {

namespace json_detail {
  template<typename>
  struct always_false : std::false_type {};

  // very minimal escaping
  inline void append_escaped(std::string& out, std::string_view sv) {
    for (char c : sv) {
      switch (c) {
        case '\\': out += "\\\\"; break;
        case '"':  out += "\\\""; break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:   out.push_back(c);
      }
    }
  }
} // namespace json_detail

template<typename T, typename = void>
struct JsonWritable : std::false_type {};

// built ins
template<>
struct JsonWritable<bool> : std::true_type {
  static void write(std::string &out, bool v) { out += v ? "true" : "false"; }
};

template<>
struct JsonWritable<std::string_view> : std::true_type {
  static void write(std::string &out, std::string_view sv) {
    out.push_back('"');
    json_detail::append_escaped(out, sv);
    out.push_back('"');
  }
};

template<typename I>
struct JsonWritable<I, std::enable_if_t<std::is_integral_v<I>>> : std::true_type {
  static void write(std::string &out, I v) { out += std::to_string(v); }
};

template<typename F>
struct JsonWritable<F, std::enable_if_t<std::is_floating_point_v<F>>> : std::true_type {
  static void write(std::string &out, F v) { out += std::to_string(v); }
};

template<typename S>
struct JsonWritable<S, std::enable_if_t<std::is_convertible_v<S, std::string_view>>> : std::true_type {
  static void write(std::string &out, S s) {
    JsonWritable<std::string_view>::write(out, std::string_view{s});
  }
};

class JsonBuilder {
  enum class Ctx { Object, Array };

  struct Frame {
    Ctx ctx;
    bool is_first_element = true;
  };

  std::string buf_;
  std::vector<Frame> stack_; // stack of contexts

  void insert_separetor() {
    auto &f = stack_.back();
    if (!f.is_first_element) buf_.push_back(',');
    f.is_first_element = false;
  }

public:
  explicit JsonBuilder(size_t reserve = 128) {
    buf_  .reserve(reserve);
    buf_  .push_back('{');
    stack_.push_back({Ctx::Object, true});
  }

  // key/value for current object
  JsonBuilder &key(std::string_view name) {
    if (stack_.empty() || stack_.back().ctx != Ctx::Object) throw std::logic_error("JsonBuilder::key() outside object");
    insert_separetor();
    JsonWritable<std::string_view>::write(buf_, name);
    buf_.push_back(':');
    return *this;
  }

  template<typename T>
  JsonBuilder &value(const T &v) {
    if (stack_.back().ctx == Ctx::Array) insert_separetor();
    JsonWritable<T>::write(buf_, v);
    return *this;
  }

  // nested structures
  JsonBuilder &begin_object() {
    if (stack_.back().ctx == Ctx::Array) insert_separetor();
    buf_  .push_back('{');
    stack_.push_back({Ctx::Object, true});
    return *this;
  }
  JsonBuilder &end_object() {
    if (stack_.empty() || stack_.back().ctx != Ctx::Object) throw std::logic_error("unbalanced end_object()");
    buf_  .push_back('}');
    stack_.pop_back();
    return *this;
  }

  JsonBuilder &begin_array() {
    if (stack_.back().ctx == Ctx::Array) insert_separetor();
    buf_  .push_back('[');
    stack_.push_back({Ctx::Array, true});
    return *this;
  }
  JsonBuilder &end_array() {
    if (stack_.empty() || stack_.back().ctx != Ctx::Array) throw std::logic_error("unbalanced end_array()");
    buf_  .push_back(']');
    stack_.pop_back();
    return *this;
  }

  // finish
  std::string str() {
    if (stack_.size() != 1 || stack_.back().ctx != Ctx::Object) throw std::logic_error("JsonBuilder::str(): JSON not closed");
    buf_  .push_back('}');
    stack_.pop_back();
    return std::move(buf_);
  }
};

class JsonReader {
  sajson::document doc_;

public:
  explicit JsonReader(std::string_view text)
      : doc_(sajson::parse(
          sajson::dynamic_allocation(),
          sajson::mutable_string_view(text.size(), const_cast<char *>(text.data())))) {
    if (!doc_.is_valid()) throw std::runtime_error(doc_.get_error_message_as_cstring());
    if (doc_.get_root().get_type() != sajson::TYPE_OBJECT) throw std::runtime_error("expected JSON object");
  }

  template<typename T> T get(std::string_view key) const {
    const sajson::string k{key.data(), key.size()};
    const sajson::value v = doc_.get_root().get_value_of_key(k);

    if (v.get_type() == sajson::TYPE_NULL)
      throw std::runtime_error("missing key \"" + std::string(key) + '"');

    return read_value<T>(v, key);
  }

  template<typename T>
  T get_or(std::string_view key, T default_value) const {
    const sajson::string k{key.data(), key.size()};
    const sajson::value v = doc_.get_root().get_value_of_key(k);

    if (v.get_type() == sajson::TYPE_NULL) return default_value;

    return read_value<T>(v, key);
  }

private:
  template<typename T>
  static T read_value(const sajson::value &v, std::string_view key) {
    if constexpr (std::is_same_v<T, bool>) {
      if (v.get_type() == sajson::TYPE_TRUE)  return true;
      if (v.get_type() == sajson::TYPE_FALSE) return false;
      throw std::runtime_error("expected boolean for key \"" + std::string(key) + '"');
    }

    else if constexpr (std::is_same_v<T, std::string>) {
      if (v.get_type() == sajson::TYPE_STRING) return v.as_string();
      throw std::runtime_error("expected string for key \"" + std::string(key) + '"');
    }

    else if constexpr (std::is_integral_v<T>) {
      if (v.get_type() == sajson::TYPE_INTEGER) return static_cast<T>(v.get_integer_value());
      if (v.get_type() == sajson::TYPE_DOUBLE)  return static_cast<T>(v.get_double_value());
      throw std::runtime_error("expected number for key \"" + std::string(key) + '"');
    }

    else if constexpr (std::is_floating_point_v<T>) {
      if (v.get_type() == sajson::TYPE_INTEGER || v.get_type() == sajson::TYPE_DOUBLE)
        return static_cast<T>(v.get_number_value());
      throw std::runtime_error("expected number for key \"" + std::string(key) + '"');
    }

    else {
      static_assert(json_detail::always_false<T>::value, "unsupported JSON type in JsonReader::get");
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_JSON_HPP
