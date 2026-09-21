#!/usr/bin/env bash
# Build, save and publish the Exorcise image for amd64 and arm64.
#
# Both architectures are published under one tag as a manifest list, so a user's
# `docker pull simonlammmm/exorcise:3.1.0` fetches whichever matches their
# machine without them naming it.
#
# Run from the repository root, not from docker/:
#
#   docker/build-and-push.sh native amd64     # on an x86 machine
#   docker/build-and-push.sh native arm64     # on an Apple silicon / Graviton machine
#   docker/build-and-push.sh manifest         # anywhere, once both are pushed
#
# Or, if you only have one machine and can wait:
#
#   docker/build-and-push.sh emulated
#
# `emulated` builds the foreign architecture under QEMU. For this image that
# means emulating a conda solve plus TensorFlow, which takes hours and
# occasionally fails outright on numeric packages. Prefer `native` where you can.

set -euo pipefail

REPO="${EXORCISE_REPO:-simonlammmm/exorcise}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
DOCKERFILE="docker/Dockerfile"

# One version for the whole project, in one file.
VERSION="$(head -n 1 "${ROOT}/VERSION" 2>/dev/null | tr -d '[:space:]')"
if [[ -z "${VERSION}" ]]; then
  echo "ERROR: could not read a version from ${ROOT}/VERSION" >&2
  exit 1
fi

# Where `save` writes tarballs, matching the existing convention in
# versions/docker/images/. This assumes the source tree sits at
# versions/src/exorcise-<version>/; set EXORCISE_IMAGE_DIR if it does not.
IMAGE_DIR="${EXORCISE_IMAGE_DIR:-$(cd "${ROOT}/../.." 2>/dev/null && pwd)/docker/images}"

usage() {
  cat <<EOF
Usage: docker/build-and-push.sh COMMAND [ARCH]

Commands:
  native ARCH   Build natively for ARCH (amd64|arm64), tag ${REPO}:${VERSION}_ARCH,
                and push it.
  save ARCH     Save an already-built ${REPO}:${VERSION}_ARCH to
                ${IMAGE_DIR}/exorcise-${VERSION}_ARCH.tar.gz
  manifest      Combine the pushed per-arch tags into ${REPO}:${VERSION}
                and ${REPO}:latest. Seconds, no rebuild.
  emulated      Build both architectures here and push with the manifest list,
                using QEMU for the foreign one. Slow.
  inspect       Show the published manifest list.

Version ${VERSION}, read from the VERSION file at the repository root.
Override the repository with EXORCISE_REPO.
EOF
}

require_arch() {
  case "${1:-}" in
    amd64 | arm64) ;;
    *) echo "ERROR: architecture must be amd64 or arm64, got '${1:-}'" >&2; exit 2 ;;
  esac
}

# A native build only needs the default builder; buildx multi-platform needs the
# docker-container driver, which the default one does not provide.
ensure_buildx_builder() {
  if ! docker buildx inspect exorcise >/dev/null 2>&1; then
    echo "Creating a buildx builder with the docker-container driver..."
    docker buildx create --name exorcise --driver docker-container
  fi
}

cmd_native() {
  require_arch "${1:-}"
  local arch="$1"
  local host
  host="$(docker version --format '{{.Server.Arch}}')"

  if [[ "${host}" != "${arch}" ]]; then
    echo "WARNING: this machine is ${host}, so building ${arch} here means QEMU" >&2
    echo "         emulation. Expect hours. Use 'emulated' if that is what you want." >&2
  fi

  cd "${ROOT}"
  docker buildx build \
    --platform "linux/${arch}" \
    -f "${DOCKERFILE}" \
    -t "${REPO}:${VERSION}_${arch}" \
    --push \
    .
  echo "Pushed ${REPO}:${VERSION}_${arch}"
}

cmd_save() {
  require_arch "${1:-}"
  local arch="$1"
  mkdir -p "${IMAGE_DIR}"
  local out="${IMAGE_DIR}/exorcise-${VERSION}_${arch}.tar.gz"

  # `docker save` needs the image locally, and a --push build does not leave it
  # there, so pull it back first.
  docker pull --platform "linux/${arch}" "${REPO}:${VERSION}_${arch}"
  docker save "${REPO}:${VERSION}_${arch}" | gzip > "${out}"
  echo "Saved ${out}"
}

cmd_manifest() {
  echo "Combining per-architecture tags into ${REPO}:${VERSION} and :latest..."
  docker buildx imagetools create \
    -t "${REPO}:${VERSION}" \
    -t "${REPO}:latest" \
    "${REPO}:${VERSION}_amd64" \
    "${REPO}:${VERSION}_arm64"
  cmd_inspect
}

cmd_emulated() {
  ensure_buildx_builder
  cd "${ROOT}"
  docker buildx build \
    --builder exorcise \
    --platform linux/amd64,linux/arm64 \
    -f "${DOCKERFILE}" \
    -t "${REPO}:${VERSION}" \
    -t "${REPO}:latest" \
    --push \
    .
  cmd_inspect
}

cmd_inspect() {
  echo
  echo "Published manifest list for ${REPO}:${VERSION}:"
  docker buildx imagetools inspect "${REPO}:${VERSION}"
}

case "${1:-}" in
  native)   shift; cmd_native "$@" ;;
  save)     shift; cmd_save "$@" ;;
  manifest) cmd_manifest ;;
  emulated) cmd_emulated ;;
  inspect)  cmd_inspect ;;
  -h | --help | help | "") usage ;;
  *) echo "Unknown command: $1" >&2; echo >&2; usage >&2; exit 2 ;;
esac
