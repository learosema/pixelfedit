/**
 * Eleventy configuration 
 * @param {import("@11ty/eleventy").UserConfig} eleventyConfig
 * @returns {import("@11ty/eleventy").UserConfig}
 */
export default function (eleventyConfig) {
  eleventyConfig.addPassthroughCopy('src/css')
  eleventyConfig.addPassthroughCopy('src/js')
  eleventyConfig.addPassthroughCopy('src/fonts')
	return {
    dir: {
      input: 'src',
      output: 'dist',
      layouts: '_layouts',
      includes: '_includes',
    }
  }
};
