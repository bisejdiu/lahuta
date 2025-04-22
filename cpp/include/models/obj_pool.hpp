#ifndef LAHUTA_OBJECT_POOL_HPP
#define LAHUTA_OBJECT_POOL_HPP

#include <cstddef>
#include <utility>
#include <vector>

// clang-format off
namespace lahuta {

template <typename T>
struct DefaultPoolTraits {
  /// creates a new instance of T (if we need to grow the pool)
  template <typename... Args>
  static T* create(Args&&... args) {
    return new T(std::forward<Args>(args)...);
  }

  /// Reset an existing T to a default/clean state
  static void reset(T* obj) {}
};

template <typename T, typename PoolTraits = DefaultPoolTraits<T>>
class ObjectPool {
public:
  explicit ObjectPool(std::size_t initial_capacity = 2000) {
    objects_.reserve(initial_capacity);
  }

  ~ObjectPool() { clear(); }

  /// Create a new object of type T, or reuse an existing one.
  inline __attribute__((always_inline))
  T* create() {
    if (next_available_index_ >= objects_.size()) {
      objects_.push_back(PoolTraits::create()); // need to grow
    }
    T* obj = objects_[next_available_index_++];
    PoolTraits::reset(obj);
    return obj;
  }

  template <typename... Args>
  T* create(Args&&... args) {
    if (next_available_index_ >= objects_.size()) {
      objects_.push_back(PoolTraits::create(std::forward<Args>(args)...));
      return objects_.back();
    }
    T* obj = objects_[next_available_index_++];
    PoolTraits::reset(obj); 
    return obj;
  }

  /// This ensures capacity for `count` objects, returning them in a vector.
  std::vector<T*> prepare(std::size_t count) {
    handle_capacity(count);
    std::vector<T*> output;
    output.reserve(count);

    next_available_index_ = 0; 
    for (std::size_t i = 0; i < count; ++i) {
      T* obj = objects_[next_available_index_++];
      PoolTraits::reset(obj);
      output.push_back(obj);
    }
    return output;
  }

  void reset() {
    next_available_index_ = 0;
  }

  // Deallocate all T objects
  void clear() {
    for (auto* ptr : objects_) {
      delete ptr;
    }
    objects_.clear();
    next_available_index_ = 0;
  }

  void handle_capacity(std::size_t required_size) {
    if (required_size > objects_.size()) {
      objects_.reserve(required_size);
      while (objects_.size() < required_size) {
        objects_.push_back(PoolTraits::create());
      }
    }
  }

private:
  std::vector<T*> objects_;
  std::size_t next_available_index_ = 0;
};

} // namespace lahuta

#endif // LAHUTA_OBJECT_POOL_HPP
