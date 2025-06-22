import { EventEmitter } from 'node:events'

/**
 * Provides some minimal APIs for the mess we do with canvas :)
 */


export class ImageData {
  /**
   * Constructor for ImageData
   * @param {Uint8ClampedArray} data
   * @param {number} width
   * @param {number} height
   */
  constructor(data, width, height) {
    this.data = data
    this.width = width
    this.height = height
  }

}

/**
 *
 * @param {Uint8ClampedArray} data
 * @param {number} width
 * @param {number} height
 * @returns
 */
export function createImageData(data, width, height) {
  return new ImageData(data, width, height)
}

/**
 * Minimal mock for Image :D
 */
export class Image {
  onload = null
  src = null
  crossOrigin = ''

  #eventEmitter = new EventEmitter()

  constructor() {
    this.#eventEmitter.on('loadimage', () => {
      if (typeof this.onload === 'function') {
        this.onload.apply(this, [this])
      }
    })
  }

  createLoadEvent() {

  }
}

export class Canvas2DRenderingContext {

}
