/**
 * Unpack bitmask data into uncompressed pixel maps (no rgba, color indexes)
 * @param {Uint8Array} bitmaskData
 * @param {number} numChars
 * @param {number} width
 * @param {number} height
 * @param {number} numBPP
 * @returns {Array<{id: number, char: string, pixels: Uint8Array}>}
 */
export function unpackFont(bitmaskData, numChars, width, height, numBPP = 1) {
  const result = []
  const bpp = Math.min(8, numBPP); // max 8 bit per pixel for now :)
  const bytesPerRow = Math.floor(0.5 + (width * bpp / 8));
  const bytesPerChar = bytesPerRow * height;
  for (let charIndex = 0; charIndex < numChars; charIndex++) {
    const pixels = new Uint8Array(width * height)
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        let val = 0;
        const xIndex = Math.floor((x * bpp) / 8);
        const bitIndex = 7 - ((x * bpp) % 8);
        const mask = ((1<<bpp)-1);
        const pixel = (bitmaskData[bytesPerChar * charIndex + bytesPerRow * y + xIndex] >> bitIndex) & mask;
        pixels[y*width + x] = pixel
      }
    }
    result.push({
      idx: charIndex,
      char: String.fromCodePoint(charIndex),
      pixels
    })
  }
  return result;
}

function packFont() {

}
