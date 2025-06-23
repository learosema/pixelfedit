import { EventEmitter } from 'node:events';

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
    this.data = data;
    this.width = width;
    this.height = height;
  }
}

/**
 * A small helper function like it is also present in the 'canvas' module
 * @param {Uint8ClampedArray} data
 * @param {number} width
 * @param {number} height
 * @returns
 */
export function createImageData(data, width, height) {
  return new ImageData(data, width, height);
}

/**
 * Minimal mock for Image :D
 */
export class Image {
  onload = null;
  src = null;
  crossOrigin = '';

  static #eventEmitter = new EventEmitter();
  shouldResolve = false


  constructor(emit = null) {
    Image.#eventEmitter.on('load', () => {
      if (typeof this.onload === 'function') {
        this.onload.apply(this, [this]);
      }
    });

    Image.#eventEmitter.on('error', (...args) => {
      if (typeof this.onerror === 'function') {
        this.onerror.apply(this, args);
      }
    });
  }

  static simulateLoad() {
    setTimeout(() =>
      Image.#eventEmitter.emit('load')
    , 0);
  }

  static simulateError(msg = 'Error') {
    setTimeout(() =>
      Image.#eventEmitter.emit('error', new Error(msg))
    , 0);
  }
}

export class Canvas2DRenderingContext {

  lineWidth = 1
  strokeStyle = '#ffffff'
  imageSmoothingEnabled = false

  #imageData = null

  constructor(width, height) {
    if (Number.isFinite(width * height)) {
      this.#imageData = new Uint8ClampedArray(width * height * 4);
    }
  }

  withImageData(imageData) {
    this.#imageData = imageData
    return this
  }

  clearRect() {}
  strokeRect() {}
  stroke() {}
  fillRect() {}
  drawImage() {}

  getImageData() {
    return this.#imageData
  }
}

export class Canvas {
  width = 400
  height = 300

  #context = new Canvas2DRenderingContext(this.width, this.height)

  getContext() {
    return this.#context
  }
}
