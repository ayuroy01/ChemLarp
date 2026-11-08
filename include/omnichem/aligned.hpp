#pragma once
// 64-byte aligned storage so SoA columns start on cache-line / AVX-512 boundaries
// and the compiler can emit aligned vector loads in hot loops.

#include <cstddef>
#include <limits>
#include <new>
#include <vector>

namespace omnichem {

inline constexpr std::size_t kSimdAlignment = 64;

template <class T, std::size_t Align = kSimdAlignment>
struct AlignedAllocator {
  using value_type = T;
  template <class U>
  struct rebind {
    using other = AlignedAllocator<U, Align>;
  };

  AlignedAllocator() noexcept = default;
  template <class U>
  AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

  [[nodiscard]] T* allocate(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) throw std::bad_array_new_length();
    return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t{Align}));
  }
  void deallocate(T* p, std::size_t) noexcept { ::operator delete(p, std::align_val_t{Align}); }

  template <class U>
  bool operator==(const AlignedAllocator<U, Align>&) const noexcept { return true; }
};

template <class T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

}  // namespace omnichem
