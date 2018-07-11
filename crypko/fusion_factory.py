import os
import time

from Crypto.Hash import keccak


class FusionFactory:
    allele_num = 2

    noise_num = 128
    bit_per_noise = 2
    noise_value_range = 2 ** bit_per_noise

    exclusive_attrs_num = 3
    bit_per_exclusive_attrs = 4
    total_exclusive_attrs_bit_num = exclusive_attrs_num * bit_per_exclusive_attrs * allele_num

    bin_attrs_num = 7
    bit_per_bin_attrs = 1

    attrs_num = exclusive_attrs_num + bin_attrs_num

    bit_per_attrs_swap = 2
    attrs_per_swap = 2 ** bit_per_attrs_swap

    bit_per_noise_mutation_pos = 4
    noise_per_mutation = 2 ** bit_per_noise_mutation_pos

    bit_per_attrs_mutation_pos = 4
    attrs_per_mutation = 2 ** bit_per_attrs_mutation_pos

    def __init__(self):
        pass

    def mix_code(self, noise1, attrs1, noise2, attrs2, entropy=None):
        """
        Given code of crypko 1 & 2, return a genetic combination - may have a random factor
        :param noise1: noise of matron
        :param attrs1: attrs of matron
        :param noise2: noise of sire
        :param attrs2: attrs of sire
        :param entropy: values used to ensure randomness
        :return: The code that are supposed to be passed down the derivative
        """
        # TODO: use global function in function body and change view to pure

        # naively and boldly use block timestamp for randomness
        # however, this should works fine, since it is one of the strongest sources of entropy in solidity
        # and the entropy of the mapping from noise/attrs to the value of the crypko is too high
        # to be manipulated by miner.
        # potentially we can split the concatenation below to get more random numbers

        if entropy is None:
            entropy = (
                round(time.time()),

                # block.blockhash(block.number - 1),  # TODO: remove `block.` for higher solidity version
                # block.coinbase,
                # block.difficulty,
            )

        rand = self.keccak(  # TODO: add back `abi.encodePacked(` for solidity v0.5+
            noise1,
            attrs1,
            noise2,
            attrs2,
            *entropy
        )

        # 1. recombine noise
        # merge origins' noise into derivative's noise
        new_noise, rand = self.recombine_batch(
            noise1,
            noise2,
            self.noise_num,
            self.bit_per_noise,
            0,
            rand
        )

        # since recombining the noise only used
        # 128 bits of the random number (at the end)
        # we can use the remaining 128 bits (at the front)
        #        require(rand >> self.noiseNum == 0);

        # recombine origins' exclusive attributes, dominant and recessive, into derivative's
        # Note: since dominant gene is on the left, and here we work reversed
        # so we recombine the recessive attributes first, and then the dominant
        # use exclusiveAttrsNum * alleleNum = 6 bits
        exclusive_attrs, rand = self.recombine(
            attrs1,
            attrs2,
            self.exclusive_attrs_num * self.allele_num,
            self.bit_per_exclusive_attrs,
            self.get_mask(self.bit_per_exclusive_attrs, 0),
            rand
        )

        # merge origins' binary attrs into derivative's
        # use binAttrsNum * alleleNum = 14 bits
        bin_attrs, rand = self.recombine_batch(
            attrs1,
            attrs2,
            self.bin_attrs_num * self.allele_num,
            self.bit_per_bin_attrs,
            self.total_exclusive_attrs_bit_num,
            rand
        )

        # since bin and exclusive attributes are on different bits,
        # we add them up to get the complete attributes
        #        require(binAttrs & exclusiveAttrs == 0);
        new_attrs = bin_attrs + exclusive_attrs

        # up until here, we right shift rand for 20 bits
        # remaining 108 bits for the rest

        # 2. swap alleles
        new_attrs, rand = self.swap_dominant_recessive(new_attrs, rand)

        # 3. mutate
        new_noise, rand = self.mutate_noise(new_noise, rand)  # ~4000 gas
        new_attrs, rand = self.mutate_attrs(new_attrs, rand)

        return new_noise, new_attrs

    def swap_dominant_recessive(self, derivative_attrs, rand):
        start_index = 0
        for i in range(self.exclusive_attrs_num):
            if rand % self.attrs_per_swap == 0:  # 25% possibility
                derivative_attrs = self.mutate(
                    derivative_attrs,
                    self.swap_bits(derivative_attrs, start_index, self.bit_per_exclusive_attrs),
                    self.get_mask(self.bit_per_exclusive_attrs * self.allele_num, start_index)
                )
            rand >>= self.bit_per_attrs_swap
            start_index += self.bit_per_exclusive_attrs * self.allele_num
        derivative_attrs, rand = self.swap_dominant_recessive_bin_attrs(derivative_attrs, rand, start_index)
        return derivative_attrs, rand

    def swap_dominant_recessive_bin_attrs(self, derivative_attrs, rand, position):
        attrs_mask = 0x5555555555555555555555555555555555555555555555555555555555555555
        rand_mask = (((rand >> 1) | rand) & attrs_mask) << position
        dominant_attrs = (((derivative_attrs >> 1) & ~rand_mask) | (derivative_attrs & rand_mask)) & attrs_mask
        rand_mask <<= 1
        recessive_attrs = (((derivative_attrs << 1) & ~rand_mask) | (derivative_attrs & rand_mask)) & (attrs_mask << 1)
        bin_attrs_mask = self.get_mask(self.bin_attrs_num * self.allele_num, position)
        derivative_attrs = self.mutate(derivative_attrs, (dominant_attrs | recessive_attrs) & bin_attrs_mask,
                                       bin_attrs_mask)
        rand >>= 2 * self.bin_attrs_num
        return derivative_attrs, rand

    def mutate_noise(self, derivative_noise, rand):
        # mutate every 16 noise
        for idx in range(0, self.noise_num, self.noise_per_mutation):  # 8 times
            # the first four bits indicate the position to mutate
            mutation_pos = (rand % self.noise_per_mutation + idx) * self.bit_per_noise
            rand >>= self.bit_per_noise_mutation_pos

            # the next two bits indicate the value to mutate
            mutation = rand % self.noise_value_range << mutation_pos
            rand >>= self.bit_per_noise

            # replace the gene at mutationPos with the mutation
            derivative_noise = self.mutate(derivative_noise, mutation, self.get_mask(self.bit_per_noise, mutation_pos))
        return derivative_noise, rand

    def mutate_attrs(self, derivative_attrs, rand):
        attr_to_mutate = rand % self.attrs_per_mutation

        # a mutation only happens on a dominant attribute.
        # since there are 14 attrs in total, we mutate just one of the attrs
        # the possibility for every attribute to mutate should be 1/16 for dominant and 0 for recessive
        if attr_to_mutate >= self.attrs_num:
            # no mutation if random num >= 14
            return derivative_attrs, rand >> self.bit_per_attrs_mutation_pos + self.bit_per_exclusive_attrs

        if attr_to_mutate < self.exclusive_attrs_num:
            # the position of the gene to mutate. it should be the allele on the left
            mutation_pos = self.allele_num * self.bit_per_exclusive_attrs * attr_to_mutate + \
                           self.bit_per_exclusive_attrs
            mutation_bit_num = self.bit_per_exclusive_attrs
        else:
            mutation_pos = self.total_exclusive_attrs_bit_num + self.allele_num * (
                    attr_to_mutate - self.exclusive_attrs_num) * self.bit_per_bin_attrs + self.bit_per_bin_attrs
            mutation_bit_num = self.bit_per_bin_attrs
        rand >>= self.bit_per_attrs_mutation_pos
        mutation_mask = self.get_mask(mutation_bit_num, mutation_pos)
        mutation = rand % 2 ** mutation_bit_num << mutation_pos
        derivative_attrs = self.mutate(derivative_attrs, mutation, mutation_mask)

        # TODO: comment line below. This is not necessary, just for code completeness.
        rand >>= self.bit_per_exclusive_attrs

        # we used at most 8 bits (or 5 if we mutate on binary attributes, or 4 if no mutation)
        # remaining 256 - 216 - 8 = 32 free bits!
        return derivative_attrs, rand

    # some helper function of computation to tidy up the codes

    @staticmethod
    def get_mask(size, position):
        # binary mask with `size` number of 1 and followed by `position` number of 0
        return 2 ** size - 1 << position

    @staticmethod
    def swap_bits(attrs, start_idx, bit_num):
        # swap the gene of size `bitNum` at `startIdx` of `attrs`
        attrs >>= start_idx
        recessive = attrs % (2 ** bit_num)
        dominant = (attrs >> bit_num) % (2 ** bit_num)
        return ((recessive << bit_num) + dominant) << start_idx

    @staticmethod
    def recombine(matron, sire, repetition, step, mask, rand):
        """
        Recombine code of matron and sire into the derivative according to the choice
        :param matron: matron's noise or attributes
        :param sire: sire's noise or attributes
        :param repetition: repeated times
        :param step: step for moving the mask
        :param mask: mask indicating the position of the currently relevant gene
        :param rand: a random number reading from right to left

        If a bit is 0, choose matron's gene, and if 1, choose sire's gene
        """
        derivative = 0
        for i in range(repetition):
            # 50% possibility choosing one of the origins
            derivative += (matron if rand % 2 == 0 else sire) & mask
            mask <<= step
            rand >>= 1
        return derivative, rand

    def recombine_batch(self, matron, sire, repetition, step, position, rand):
        """
        Recombine code of matron and sire into the derivative according to the choice
        :param matron: matron's noise or attributes
        :param sire: sire's noise or attributes
        :param repetition: repeated times
        :param step: step for moving the mask
        :param position: indicating the starting position
        :param rand: a random number reading from right to left

        If a bit is 0, choose matron's gene, and if 1, choose sire's gene
        """
        rand_dup = self.duplicate_bits(rand, step)
        rand_dup <<= position
        mask = self.get_mask(step * repetition, position)
        derivative = ((~rand_dup & matron) | (rand_dup & sire)) & mask
        rand >>= repetition
        return derivative, rand

    @staticmethod
    def duplicate_bits(input_val, count):
        """
        Duplicate every bit for `count` times, ignore the highest bits
        Only works when count == 1 (do nothing) or count == 2

        >>> FusionFactory.duplicate_bits(0b01101011, 2)
        11001111

        :param input_val: The input number
        :param count: The number duplicate bits for each input bit
        """
        if count == 1:
            return input_val
        elif count == 2:
            mask1 = 0x00000000000000000000000000000000ffffffffffffffff0000000000000000
            mask2 = 0x000000000000000000000000000000000000000000000000ffffffffffffffff
            input_val = ((input_val & mask1) << 64) | (input_val & mask2)
            mask1 = 0x0000000000000000ffffffff000000000000000000000000ffffffff00000000
            mask2 = 0x000000000000000000000000ffffffff000000000000000000000000ffffffff
            input_val = ((input_val & mask1) << 32) | (input_val & mask2)
            mask1 = 0x00000000ffff000000000000ffff000000000000ffff000000000000ffff0000
            mask2 = 0x000000000000ffff000000000000ffff000000000000ffff000000000000ffff
            input_val = ((input_val & mask1) << 16) | (input_val & mask2)
            mask1 = 0x0000ff000000ff000000ff000000ff000000ff000000ff000000ff000000ff00
            mask2 = 0x000000ff000000ff000000ff000000ff000000ff000000ff000000ff000000ff
            input_val = ((input_val & mask1) << 8) | (input_val & mask2)
            mask1 = 0x00f000f000f000f000f000f000f000f000f000f000f000f000f000f000f000f0
            mask2 = 0x000f000f000f000f000f000f000f000f000f000f000f000f000f000f000f000f
            input_val = ((input_val & mask1) << 4) | (input_val & mask2)
            mask1 = 0x0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c
            mask2 = 0x0303030303030303030303030303030303030303030303030303030303030303
            input_val = ((input_val & mask1) << 2) | (input_val & mask2)
            mask1 = 0x2222222222222222222222222222222222222222222222222222222222222222
            mask2 = 0x1111111111111111111111111111111111111111111111111111111111111111
            input_val = ((input_val & mask1) << 1) | (input_val & mask2)
            input_val = (input_val << 1) | input_val
            return input_val
        else:
            raise ValueError  # not supported

    @staticmethod
    def mutate(derivative, mutation, mask):
        """
        Mutate the gene at the `mask` position of the derivative code to `mutation`
        """
        return derivative & ~mask + mutation

    @staticmethod
    def keccak(*args):
        keccak_hash = keccak.new(digest_bits=256)
        for i in args:
            keccak_hash.update(i.to_bytes(len(bin(i)) - 1, 'big'))
        return int(keccak_hash.hexdigest(), 16)


def binary_average(values, bits):
    bitmap = [
        [(j & (2 ** i)) >> i for j in values] for i in range(bits)
    ]
    averages = [
        round(sum(i) / len(i)) for i in bitmap
    ]
    return sum(
        i << n for n, i in enumerate(averages)
    )


def approximate(noise1, attrs1, noise2, attrs2, average=1):
    attrs = []

    factory = FusionFactory()
    for _ in range(average):
        attrs.append(factory.mix_code(noise1, attrs1, noise2, attrs2,
                                      entropy=(int.from_bytes(os.urandom(64), 'big'),))[1])

    return binary_average(attrs, 38)


if __name__ == '__main__':
    n1 = 90473127252118133248908235787953752405933676411820653807769447369544567230866
    a1 = 52818018743
    n2 = 3741807294364386106505680355724161895371676705092164253676086926199092543794
    a2 = 52818020023

    print(f'In [0]: {bin(a1)}')
    print(f'In [1]: {bin(a2)}')

    print()
    for i in range(10):
        na3 = approximate(n1, a1, n2, a2)
        print(f'Out[{i}]: {bin(na3)}')

    n3 = 90587641027505316411294546399164831736118732119729757568075160231330036711826
    a3 = 52818022407

    print(f'\nOut[0]: {bin(a3)}')
