import { describe, it } from 'node:test'
import assert from 'node:assert/strict'

import { unpackFont } from "./packer.js"

describe('packer', () => {
  describe('unpackFont', () => {

    it('should unpack basic bitmask data into pixels', () => {
      const font = unpackFont(new Uint8Array([
        0b00000000,
        0b01111100,
        0b01100110,
        0b01111100,
        0b01100110,
        0b01100110,
        0b01111100,
        0b00000000,
      ]), 1, 8, 8, 1)
      assert.equal(font.length, 1)
      assert.deepEqual(font[0].pixels, new Uint8Array([
        0,0,0,0,0,0,0,0,
        0,1,1,1,1,1,0,0,
        0,1,1,0,0,1,1,0,
        0,1,1,1,1,1,0,0,
        0,1,1,0,0,1,1,0,
        0,1,1,0,0,1,1,0,
        0,1,1,1,1,1,0,0,
        0,0,0,0,0,0,0,0,
      ]))
    })

  })
})
