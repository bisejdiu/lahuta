#ifndef LAHUTA_PROPS_HPP
#define LAHUTA_PROPS_HPP

namespace lahuta {

struct LuniProperties {
  static void register_all();

  static void initialize() {
    static bool initialized = (register_all(), true);
    (void)initialized;
  }
};

} // namespace lahuta

#endif // LAHUTA_PROPS_HPP
