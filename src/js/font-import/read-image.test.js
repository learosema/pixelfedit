import { describe, it } from 'node:test';
import assert, { deepStrictEqual } from 'node:assert';

import { createImageData } from '../testing/imagemocks.js';
import { convertImageData, getImageColors } from './read-image.js';

function generateTestPattern(numColors, width, height) {
  return createImageData(
    new Uint8ClampedArray(Array.from({length: width * height}, (_, idx) => [idx % numColors,0,0,255]).flat()),
    width, height
  )
}

function createTestImageData(str) {
  const lines = str.trim().split('\n').map(str => str.trim().replace(/\s/g, ''))
  const width = Math.max(...lines.map(l => l.length))
  const height = lines.length;
  const data =  new Uint8ClampedArray(Array.from({length: width * height}))
  for (let y = 0; y < height; y++) {
    const line = lines[y]
    if (! line) {
      continue;
    }
    for (let x = 0; x < width; x++) {
      const color = parseInt(line[x], 10) || 0
      data[(y * width + x) * 4] = color
      data[(y * width + x) * 4 + 1] = 0
      data[(y * width + x) * 4 + 2] = 0
      data[(y * width + x) * 4 + 3] = 255
    }
  }
  return createImageData(data, width, height);
}



describe('getImageColors', () => {

  it('should count colors correctly', () => {
    const testData = generateTestPattern(5, 16, 16)
    const colors = getImageColors(testData)
    assert.deepStrictEqual(colors, ["#000000ff", "#000001ff", "#000002ff", "#000003ff", "#000004ff"])
  })

  it('should count colors correctly but limit it if specified', () => {
    const testData = generateTestPattern(5, 16, 16)
    const colors = getImageColors(testData, 2)
    assert.deepStrictEqual(colors, ["#000000ff", "#000001ff"])
  })
})


describe('convertImageData', () => {
  const testData = createTestImageData(`
    ........ ........
    .11111.. ..1111..
    .1....1. .1....1.
    .1....1. .1....1.
    .11111.. .1...11.
    .1...... ..11111.
    .1...... .....111
    ........ ........
  `)

  it('should convert image data to bitmask data', () => {

    const result = convertImageData(testData, 8, 8)

    const P = [
      0b00000000,
      0b01111100,
      0b01000010,
      0b01000010,
      0b01111100,
      0b01000000,
      0b01000000,
      0b00000000
    ]

    const Q = [
      0b00000000,
      0b00111100,
      0b01000010,
      0b01000010,
      0b01000110,
      0b00111110,
      0b00000011,
      0b00000000
    ]

    deepStrictEqual(result, [...P, ...Q])
  })
})

